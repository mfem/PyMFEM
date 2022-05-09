'''
   MFEM example 6
      This is a version of Example 1 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
from mpi4py import MPI
import numpy as np
import os

num_proc = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)
verbose = (myid == 0)


def run(meshfile="",
        order=1,
        smooth_el=True,
        max_dofs=100000,
        reorder_mesh=0,
        restart=False,
        nc_simplices=True,
        visualization=False,
        device_config="cpu",
        pa=False):

    device = mfem.Device(device_config)
    if myid == 0:
        device.Print()

    # 4. Read the (serial) mesh from the given mesh file on all processors.
    #    We can handle triangular, quadrilateral, tetrahedral, hexahedral,
    #    surface and volume meshes with the same code.

    mesh = mfem.Mesh(meshfile, 1, 1)

    # 5. A NURBS mesh cannot be refined locally so we refine it uniformly
    #    and project it to a standard curvilinear mesh of order 2.
    if (mesh.NURBSext):
        mesh.UniformRefinement()
        mesh.SetCurvature(2)

    # 6. MFEM supports dynamic partitioning (load balancing) of parallel non-
    #    conforming meshes based on space-filling curve (SFC) partitioning.
    #    SFC partitioning is extremely fast and scales to hundreds of
    #    thousands of processors, but requires the coarse mesh to be ordered,
    #    ideally as a sequence of face-neighbors. The mesh may already be
    #    ordered (like star-hilbert.mesh) or we can order it here. Ordering
    #    type 1 is a fast spatial sort of the mesh, type 2 is a high quality
    #    optimization algorithm suitable for ordering general unstructured
    #    meshes.
    print("reorder_mesh", reorder_mesh)
    if reorder_mesh > 0:
        ordering = mfem.intArray()
        if reorder_mesh == 1:
            mesh.GetHilbertElementOrdering(ordering)
        elif reorder_mesh == 2:
            mesh.GetGeckoElementOrdering(ordering)
        mesh.ReorderElements(ordering)

    # 7. Make sure the mesh is in the non-conforming mode to enable local
    #    refinement of quadrilaterals/hexahedra, and the above partitioning
    #    algorithm. Simplices can be refined either in conforming or in non-
    #    conforming mode. The conforming mode however does not support
    #    dynamic partitioning.
    mesh.EnsureNCMesh(nc_simplices)

    # 8. Define a parallel mesh by partitioning the serial mesh.
    #    Once the parallel mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)

    dim = pmesh.Dimension()
    sdim = pmesh.SpaceDimension()

    ess_bdr = intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)

    # 10. Define a finite element space on the mesh. The polynomial order is
    #     one (linear) by default, but this can be changed on the command line.
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)

    # 11. As in Example 1p, we set up bilinear and linear forms corresponding to
    #     the Laplace problem -\Delta u = 1. We don't assemble the discrete
    #     problem yet, this will be done in the main loop.
    a = mfem.ParBilinearForm(fespace)
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
        a.SetDiagonalPolicy(mfem.Operator.DIAG_ONE)
    b = mfem.ParLinearForm(fespace)

    one = mfem.ConstantCoefficient(1.0)

    integ = mfem.DiffusionIntegrator(one)
    a.AddDomainIntegrator(integ)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))

    # 12. The solution vector x and the associated finite element grid function
    #     will be maintained over the AMR iterations. We initialize it to zero.
    x = mfem.ParGridFunction(fespace)
    x.Assign(0.0)

    # 13. Connect to GLVis.
    if visualization:
        sout = mfem.socketstream('localhost', 19916)
        if sout.good():
            sout.precision(8)
        else:
            visualization = False

    # 14. Set up an error estimator. Here we use the Zienkiewicz-Zhu estimator
    #     with L2 projection in the smoothing step to better handle hanging
    #     nodes and parallel partitioning. We need to supply a space for the
    #     discontinuous flux (L2) and a space for the smoothed flux

    flux_fec = mfem.L2_FECollection(order, dim)
    flux_fes = mfem.ParFiniteElementSpace(pmesh, flux_fec, sdim)

    if smooth_el and dim > 1:
        smooth_flux_fec = mfem.RT_FECollection(order-1, dim)
        smooth_flux_fes = mfem.ParFiniteElementSpace(pmesh, smooth_flux_fec, 1)
    else:
        smooth_flux_fec = mfem.H1_FECollection(order, dim)
        smooth_flux_fes = mfem.ParFiniteElementSpace(
            pmesh, smooth_flux_fec, dim)

    estimator = mfem.L2ZienkiewiczZhuEstimator(integ, x, flux_fes,
                                               smooth_flux_fes)
    # 15. A refiner selects and refines elements based on a refinement strategy.
    #     The strategy here is to refine elements with errors larger than a
    #     fraction of the maximum element error. Other strategies are possible.
    #     The refiner will call the given error estimator.
    refiner = mfem.ThresholdRefiner(estimator)
    refiner.SetTotalErrorFraction(0.7)

    # 16. The main AMR loop. In each iteration we solve the problem on the
    #     current mesh, visualize the solution, and refine the mesh.
    it = 0
    while True:
        global_dofs = fespace.GlobalTrueVSize()
        if (myid == 0):
            print("AMR iteration " + str(it))
            print("Number of unknowns: " + str(global_dofs))

        # 17. Assemble the right-hand side and determine the list of true
        #     (i.e. parallel conforming) essential boundary dofs.
        ess_tdof_list = intArray()
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)
        b.Assemble()

        # 18. Assemble the stiffness matrix. Note that MFEM doesn't care at this
        #     point that the mesh is nonconforming and parallel.  The FE space is
        #     considered 'cut' along hanging edges/faces, and also across
        #     processor boundaries.
        a.Assemble()

        # 19. Create the linear system: eliminate boundary conditions.
        #     The system will be solved for true (unconstrained) DOFs only.
        A = mfem.OperatorPtr()
        B = mfem.Vector()
        X = mfem.Vector()
        copy_interior = 1
        a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior)

        # 20. Solve the linear system A X = B.
        #     * With full assembly, use the BoomerAMG preconditioner from hypre.
        #     * With partial assembly, use a diagonal preconditioner.
        if pa:
            M = mfem.OperatorJacobiSmoother(a, ess_tdof_list)
        else:
            amg = mfem.HypreBoomerAMG()
            amg.SetPrintLevel(0)
            M = amg

        pcg = mfem.CGSolver(MPI.COMM_WORLD)
        pcg.SetPreconditioner(M)
        pcg.SetOperator(A.Ptr())
        pcg.SetRelTol(1e-6)
        pcg.SetMaxIter(2000)
        pcg.SetPrintLevel(3)
        pcg.Mult(B, X)

        # 21. Switch back to the host and extract the parallel grid function
        #     corresponding to the finite element approximation X. This is the
        #     local solution on each processor.
        a.RecoverFEMSolution(X, b, x)
        if (global_dofs > max_dofs):
            if (myid == 0):
                print("Reached the maximum number of dofs. Stop.")
            break

        # 22. Send the solution by socket to a GLVis server.
        if visualization:
            sout << "parallel " << num_proc << " " << myid << "\n"
            sout << "solution\n" << pmesh << x
            sout.flush()

        # 23. Call the refiner to modify the mesh. The refiner calls the error
        #     estimator to obtain element errors, then it selects elements to be
        #     refined and finally it modifies the mesh. The Stop() method can be
        #     used to determine if a stopping criterion was met.

        refiner.Apply(pmesh)
        if (refiner.Stop()):
            if myid == 0:
                print("Stopping criterion satisfied. Stop.")
            break
        # 24. Update the finite element space (recalculate the number of DOFs,
        #     etc.) and create a grid function update matrix. Apply the matrix
        #     to any GridFunctions over the space. In this case, the update
        #     matrix is an interpolation matrix so the updated GridFunction will
        #     still represent the same function as before refinement.
        fespace.Update()
        x.Update()

        # 25. Load balance the mesh, and update the space and solution. Currently
        #     available only for nonconforming meshes.
        if pmesh.Nonconforming():
            pmesh.Rebalance()
            # Update the space and the GridFunction. This time the update matrix
            # redistributes the GridFunction among the processors.
            fespace.Update()
            x.Update()

        # 26. Inform also the bilinear and linear forms that the space has
        #     changed.
        a.Update()
        b.Update()

        # 27. Save the current state of the mesh every 5 iterations. The
        #     computation can be restarted from this point. Note that unlike in
        #     visualization, we need to use the 'ParPrint' method to save all
        #     internal parallel data structures.

        if ((it + 1) % 5) == 0:
            pmesh.ParPrint("ex6p-checkpoint."+smyid, 8)
            if myid == 0:
                print("Checkpoint saved")
        it = it + 1


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex6p (Laplace problem with AMR)')
    parser.add_argument('-m', '--mesh',
                        default="star.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store',
                        default=1,
                        type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    parser.add_argument("-rm", "--reorder-mesh",
                        action='store',
                        default=0,
                        type=int,
                        help="".join(["Reorder elements of the coarse mesh to improve ",
                                      "dynamic partitioning: 0=none, 1=hilbert, 2=gecko."]))
    parser.add_argument("-cs", "--conforming-simplices",
                        action='store_true',
                        default=False,
                        help="For simplicial meshes, disable nonconforming refinement")
    parser.add_argument("-md", "--max-dofs",
                        action="store",
                        type=int,
                        default=100000,
                        help="Stop after reaching this many degrees of freedom.")
    parser.add_argument("-h1", "--smooth-h1",
                        action='store_true',
                        default=False,
                        help="Represent the smooth flux in vector H1 space.")
    # parser.add_argument("-res", "--restart",
    #                    action='store_true',
    #                    help="Restart computation from the last checkpoint.");

    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    smooth_el = not args.smooth_h1

    run(meshfile=meshfile,
        order=args.order,
        smooth_el=smooth_el,
        max_dofs=args.max_dofs,
        reorder_mesh=args.reorder_mesh,
        restart=False,
        nc_simplices=not args.conforming_simplices,
        visualization=not args.no_visualization,
        device_config=args.device,
        pa=args.partial_assembly)
