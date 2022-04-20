'''
   MFEM example 21
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(order=1,
        static_cond=False,
        meshfile="",
        visualization=1,
        serial_ref_levels=0):

    device = mfem.Device('cpu')
    if myid == 0:    
        device.Print()
    
    # 2. Read the mesh from the given mesh file. We can handle triangular,
    #    quadrilateral, tetrahedral, and hexahedral meshes with the same code.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    assert mesh.SpaceDimension() == dim, "invalid mesh"

    if mesh.attributes.Max() < 2 or mesh.bdr_attributes.Max() < 2:
        print("\n".join(["Input mesh should have at least two materials and "
                         "two boundary attributes! (See schematic in ex2.cpp)"]))

    # 3. Since a NURBS mesh can currently only be refined uniformly, we need to
    #    convert it to a piecewise-polynomial curved mesh. First we refine the
    #    NURBS mesh a bit more and then project the curvature to quadratic Nodes.
    if mesh.NURBSext and serial_ref_levels == 0:
        serial_ref_levels = 2

    for i in range(serial_ref_levels):
        mesh.UniformRefinement()

    if mesh.NURBSext:
        mesh.SetCurvature(2)
    mesh.EnsureNCMesh()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)

    # 4. Define a finite element space on the mesh. The polynomial order is
    #    one (linear) by default, but this can be changed on the command line.
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec, dim)

    # 5. As in Example 2, we set up the linear form b(.) which corresponds to
    #    the right-hand side of the FEM linear system. In this case, b_i equals
    #    the boundary integral of f*phi_i where f represents a "pull down"
    #    force on the Neumann part of the boundary and phi_i are the basis
    #    functions in the finite element fespace. The force is defined by the
    #    VectorArrayCoefficient object f, which is a vector of Coefficient
    #    objects. The fact that f is non-zero on boundary attribute 2 is
    #    indicated by the use of piece-wise constants coefficient for its last
    #    component. We don't assemble the discrete problem yet, this will be
    #    done in the main loop.
    f = mfem.VectorArrayCoefficient(dim)
    for i in range(dim-1):
        f.Set(i, mfem.ConstantCoefficient(0.0))

    pull_force = mfem.Vector(pmesh.bdr_attributes.Max())
    pull_force.Assign(0.0)
    pull_force[1] = -1.0e-2
    f.Set(dim-1, mfem.PWConstCoefficient(pull_force))

    b = mfem.ParLinearForm(fespace)
    b.AddDomainIntegrator(mfem.VectorBoundaryLFIntegrator(f))

    # 6. Set up the bilinear form a(.,.) on the finite element space
    #    corresponding to the linear elasticity integrator with piece-wise
    #    constants coefficient lambda and mu.
    llambda = mfem.Vector(pmesh.attributes.Max())
    llambda.Assign(1.0)
    llambda[0] = llambda[1]*50
    lambda_func = mfem.PWConstCoefficient(llambda)
    mu = mfem.Vector(pmesh.attributes.Max())
    mu.Assign(1.0)
    mu[0] = mu[1]*50
    mu_func = mfem.PWConstCoefficient(mu)

    a = mfem.ParBilinearForm(fespace)
    integ = mfem.ElasticityIntegrator(lambda_func, mu_func)
    a.AddDomainIntegrator(integ)
    if static_cond:
        a.EnableStaticCondensation()

    # 7. The solution vector x and the associated finite element grid function
    #    will be maintained over the AMR iterations. We initialize it to zero.
    zero_vec = mfem.Vector(dim)
    zero_vec.Assign(0.0)
    zero_vec_coeff = mfem.VectorConstantCoefficient(zero_vec)
    x = mfem.ParGridFunction(fespace)
    x.Assign(0.0)

    # 8. Determine the list of true (i.e. conforming) essential boundary dofs.
    #    In this example, the boundary conditions are defined by marking only
    #    boundary attribute 1 from the mesh as essential and converting it to a
    #    list of true dofs.  The conversion to true dofs will be done in the
    #    main loop.
    ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(0)
    ess_bdr[0] = 1

    # 9. Connect to GLVis.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)

    # 10. Set up an error estimator. Here we use the Zienkiewicz-Zhu estimator
    #     that uses the ComputeElementFlux method of the ElasticityIntegrator to
    #     recover a smoothed flux (stress) that is subtracted from the element
    #     flux to get an error indicator. We need to supply the space for the
    #     smoothed flux: an (H1)^tdim (i.e., vector-valued) space is used here.
    #     Here, tdim represents the number of components for a symmetric (dim x
    #     dim) tensor.
    tdim = dim*(dim+1)//2
    flux_fec = mfem.L2_FECollection(order, dim)
    flux_fespace = mfem.ParFiniteElementSpace(pmesh, flux_fec, tdim)
    smooth_flux_fespace = mfem.ParFiniteElementSpace(pmesh, fec, tdim)
    estimator = mfem.L2ZienkiewiczZhuEstimator(
        integ, x, flux_fespace, smooth_flux_fespace)

    # 11. A refiner selects and refines elements based on a refinement strategy.
    #     The strategy here is to refine elements with errors larger than a
    #     fraction of the maximum element error. Other strategies are possible.
    #     The refiner will call the given error estimator.
    refiner = mfem.ThresholdRefiner(estimator)
    refiner.SetTotalErrorFraction(0.7)

    # 12. The main AMR loop. In each iteration we solve the problem on the
    #     current mesh, visualize the solution, and refine the mesh.
    max_dofs = 50000
    max_amr_itr = 20
    for it in range(max_amr_itr+1):
        global_dofs = fespace.GlobalTrueVSize()
        if myid == 0:
            print("\n".join(["",
                             "AMR iteration " + str(it),
                             "Number of unknowns: " + str(global_dofs)]))

        # 13. Assemble the stiffness matrix and the right-hand side.
        a.Assemble()
        b.Assemble()

        # 14. Set Dirichlet boundary values in the GridFunction x.
        #     Determine the list of Dirichlet true DOFs in the linear system.
        ess_tdof_list = mfem.intArray()
        x.ProjectBdrCoefficient(zero_vec_coeff, ess_bdr)
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

        # 15. Create the linear system: eliminate boundary conditions, constrain
        #     hanging nodes and possibly apply other transformations. The system
        #     will be solved for true (unconstrained) DOFs only.
        A = mfem.HypreParMatrix()
        X = mfem.Vector()
        B = mfem.Vector()
        copy_interior = 1
        a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior)

        # 16. Define and apply a parallel PCG solver for AX=B with the BoomerAMG
        #     preconditioner from hypre.
        amg = mfem.HypreBoomerAMG()
        amg.SetPrintLevel(0)
        # amg.SetSystemsOptions(dim); // optional
        pcg = mfem.CGSolver(A.GetComm())
        pcg.SetPreconditioner(amg)
        pcg.SetOperator(A)
        pcg.SetRelTol(1e-6)
        pcg.SetMaxIter(500)
        pcg.SetPrintLevel(3)  # print the first and the last iterations only
        pcg.Mult(B, X)

        # 17. After solving the linear system, reconstruct the solution as a
        #     finite element GridFunction. Constrained nodes are interpolated
        #     from true DOFs (it may therefore happen that x.Size() >= X.Size()).
        a.RecoverFEMSolution(X, b, x)

        # 18. Send solution by socket to the GLVis server.
        if visualization and sol_sock.good():
            nodes = mfem.GridFunction(fespace)
            pmesh.GetNodes(nodes)
            nodes += x
            own_nodes = 0
            nodes, own_nodes = pmesh.SwapNodes(nodes, own_nodes)

            x.Neg()  # visualize the backward displacement
            sol_sock << "parallel " << num_procs << ' ' << myid << '\n'
            sol_sock << "solution\n" << mesh << x
            sol_sock.flush()
            x.Neg()

            nodes, own_nodes = pmesh.SwapNodes(nodes, own_nodes)
            if it == 0:
                kk = "Rj1" if dim == 2 else ""
                sol_sock << "keys '" << kk << "m'" << "\n"
                sol_sock.endline()

            sol_sock << "window_title 'AMR iteration: " << it << "'\n" << "pause"

            if myid == 0:
                print("Visualization paused. "
                      "Press <space> in the GLVis window to continue.")

        if (global_dofs > max_dofs):
            if myid == 0:
                print("Reached the maximum number of dofs. Stop.")
            break

        # 19. Call the refiner to modify the mesh. The refiner calls the error
        #     estimator to obtain element errors, then it selects elements to be
        #     refined and finally it modifies the mesh. The Stop() method can be
        #     used to determine if a stopping criterion was met.
        refiner.Apply(pmesh)
        if (refiner.Stop()):
            if myid == 0:
                print("Stopping criterion satisfied. Stop.")
            break

        # 20. Update the space to reflect the new state of the mesh. Also,
        #     interpolate the solution x so that it lies in the new space but
        #     represents the same function. This saves solver iterations later
        #     since we'll have a good initial guess of x in the next step.
        #     Internally, FiniteElementSpace::Update() calculates an
        #     interpolation matrix which is then used by GridFunction::Update().
        fespace.Update()
        x.Update()

        # 21. Load balance the mesh, and update the space and solution. Currently
        #     available only for nonconforming meshes.
        if pmesh.Nonconforming():
            pmesh.Rebalance()
            fespace.Update()
            x.Update()
        a.Update()
        b.Update()

    pmesh.Print("ex21p_reference_mesh"+smyid, 16)
    nodes = mfem.GridFunction(fespace)
    pmesh.GetNodes(nodes)
    nodes += x
    own_nodes = 0
    nodes2, own_nodes = pmesh.SwapNodes(nodes, own_nodes)
    pmesh.Print("ex21p_deformed_mesh"+smyid, 16)
    nodes, own_nodes = pmesh.SwapNodes(nodes2, own_nodes)

    x.Save("ex21p_displacement"+smyid)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex21 (Adaptive mesh refinement for linear elasticity)')
    parser.add_argument('-m', '--mesh',
                        default='beam-tri.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-rs", "--refine_serial",
                        action='store', type=int, default=0,
                        help="Number of serial refinement before parallel partitioning")
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-sc', '--static-condensation',
                        action='store_true',
                        help="Enable static condensation.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=True,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    order = args.order
    static_cond = args.static_condensation

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization

    run(order=order,
        static_cond=static_cond,
        meshfile=meshfile,
        visualization=visualization,
        serial_ref_levels=args.refine_serial)
