'''
   PyMFEM example 34p

   See c++ version in the MFEM library for more detail

   Sample runs:  mpirun -np 4 python ex34p.py -o 2
                 mpirun -np 4 python ex34p.py -pa -hex

   Device sample runs:
                 mpirun -np 4 python ex34p.py -o 2 -pa -hex -d cuda
                 mpirun -np 4 python ex34p.py -o 2 -no-pa -d cuda

'''
import mfem.par as mfem
from mfem.par import intArray, doubleArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.{:0>6d}'.format(myid)

pa = False
algebraic_ceed = False
visualization = True

def run(ser_ref_levels=1,
        par_ref_levels=1,
        order=1,
        delta_const=1e-6,
        static_cond=False,
        mixed=True,
        device_config="cpu"):

    cond_attr = mfem.intArray()

    mesh_file = "../data/fichera-mixed.mesh"
    if not mixed or pa:
        mesh_file = "../data/fichera.mesh"

    if mesh_file == "../data/fichera-mixed.mesh":
        submesh_elems = mfem.intArray([0, 2, 3, 4, 9])
    else:
        submesh_elems = mfem.intArray([10, 14, 34, 36, 37, 38, 39])

    sym_plane_attr = mfem.intArray([9, 10, 11, 12, 13, 14, 15, 16, ])

    phi0_attr = mfem.intArray()
    phi1_attr = mfem.intArray()
    jn_zero_attr = mfem.intArray()

    if (mesh_file == "../data/fichera-mixed.mesh" or
            mesh_file == "../data/fichera.mesh"):
        phi0_attr.Append(2)

    if (mesh_file == "../data/fichera-mixed.mesh" or
            mesh_file == "../data/fichera.mesh"):
        phi1_attr.Append(23)

    if (mesh_file == "../data/fichera-mixed.mesh" or
            mesh_file == "../data/fichera.mesh"):
        jn_zero_attr.Append(25)
    for i in range(sym_plane_attr.Size()):
        jn_zero_attr.Append(sym_plane_attr[i])

    # 3. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device_config)
    device.Print()

    # 4. Read the (serial) mesh from the given mesh file on all processors.  We
    #    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    #    and volume meshes with the same code.
    mesh_file = expanduser(
        join(os.path.dirname(__file__), *mesh_file.split("/")))
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()
    if not mixed or pa:
        mesh.UniformRefinement()
        if ser_ref_levels > 0:
            ser_ref_levels = ser_ref_levels - 1
        else:
            par_ref_levels = par_ref_levels - 1
    submesh_attr = -1

    if cond_attr.Size() == 0 and submesh_elems.Size() > 0:
        max_attr = mesh.attributes.Max()
        submesh_attr = max_attr + 1

        for i in range(submesh_elems.Size()):
            mesh.SetAttribute(submesh_elems[i], submesh_attr)
        mesh.SetAttributes()
        if cond_attr.Size() == 0:
            cond_attr.Append(submesh_attr)

    # 5. Refine the serial mesh on all processors to increase the resolution. In
    #    this example we do 'ref_levels' of uniform refinement.

    for l in range(ser_ref_levels):
        mesh.UniformRefinement()

    # 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution. Once the
    #    parallel mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh
    for l in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 6b. Extract a submesh covering a portion of the domain
    pmesh_cond = mfem.ParSubMesh.CreateFromDomain(pmesh, cond_attr)

    # 7. Define a suitable finite element space on the SubMesh and compute
    #    the current density as an H(div) field.
    fec_cond_rt = mfem.RT_FECollection(order - 1, dim)
    fes_cond_rt = mfem.ParFiniteElementSpace(pmesh_cond, fec_cond_rt)
    j_cond = mfem.ParGridFunction(fes_cond_rt)

    ComputeCurrentDensityOnSubMesh(
        order, phi0_attr, phi1_attr, jn_zero_attr, j_cond)

    # 7a. Save the SubMesh and associated current density in parallel. This
    #     output can be viewed later using GLVis:
    #        "glvis -np <np> -m cond_mesh -g cond_j"
    pmesh_cond.Print("cond_mesh"+smyid, 8)
    j_cond.Save("cond_j"+smyid, 8)

    # 7b. Send the current density, computed on the SubMesh, to a GLVis server.
    if visualization:
        port_sock = mfem.socketstream("localhost", 19916)
        port_sock << "parallel " << num_procs << " " << myid << "\n"
        port_sock.precision(8)
        port_sock << "solution\n" << pmesh_cond << j_cond
        port_sock << "window_title 'Conductor J'"
        port_sock << "window_geometry 400 0 400 350"
        port_sock.flush()

    # 8. Define a parallel finite element space on the full mesh. Here we use
    #    the H(curl) finite elements for the vector potential and H(div) for the
    #    current density.
    fec_nd = mfem.ND_FECollection(order, dim)
    fec_rt = mfem.RT_FECollection(order - 1, dim)
    fespace_nd = mfem.ParFiniteElementSpace(pmesh, fec_nd)
    fespace_rt = mfem.ParFiniteElementSpace(pmesh, fec_rt)

    j_full = mfem.ParGridFunction(fespace_rt)
    j_full.Assign(0.0)
    pmesh_cond.Transfer(j_cond, j_full)

    # 8a. Send the transferred current density to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << pmesh << j_full
        sol_sock << "window_title 'J Full'"
        sol_sock << "window_geometry 400 430 400 350"
        sol_sock.flush()

    # 9. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, the boundary conditions are defined by
    #    marking all the boundary attributes except for those on a symmetry
    #    plane as essential (Dirichlet) and converting them to a list of true
    #    dofs.
    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray()
    if pmesh.bdr_attributes.Size():
        ess_bdr.SetSize(pmesh.bdr_attributes.Max())
        ess_bdr.Assign(1)
        for i in range(sym_plane_attr.Size()):
            ess_bdr[sym_plane_attr[i]-1] = 0
        fespace_nd.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 10. Set up the parallel linear form b(.) which corresponds to the
    #     right-hand side of the FEM linear system, which in this case is
    #     (J,W_i) where J is given by the function H(div) field transferred
    #     from the SubMesh and W_i are the basis functions in the finite
    #     element fespace.
    jCoef = mfem.VectorGridFunctionCoefficient(j_full)
    b = mfem.ParLinearForm(fespace_nd)
    b.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(jCoef))
    b.Assemble()

    # 11. Define the solution vector x as a parallel finite element grid
    #     function corresponding to fespace. Initialize x to zero.
    x = mfem.ParGridFunction(fespace_nd)
    x.Assign(0.0)

    # 12. Set up the parallel bilinear form corresponding to the EM diffusion
    #     operator curl muinv curl + delta I, by adding the curl-curl and the
    #     mass domain integrators. For standard magnetostatics equations choose
    #     delta << 1. Larger values of delta should make the linear system
    #     easier to solve at the expense of resembling a diffusive quasistatic
    #     magnetic field.  A reasonable balance must be found whenever the mesh
    #     or problem setup is altered.
    muinv = mfem.ConstantCoefficient(1.0)
    delta = mfem.ConstantCoefficient(delta_const)
    a = mfem.ParBilinearForm(fespace_nd)
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(delta))

    # 13. Assemble the parallel bilinear form and the corresponding linear
    #     system, applying any necessary transformations such as: parallel
    #     assembly, eliminating boundary conditions, applying conforming
    #     constraints for non-conforming AMR, static condensation, etc.
    if static_cond:
        a.EnableStaticCondensation()
    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

    # 14. Solve the system AX=B using PCG with an AMS preconditioner.
    if myid == 0:
        print("\n".join(("\nSolving for magnetic vector potential ",
                         "using CG with AMS")))

    if pa:  # Jacobi preconditioning in partial assembly mode
        ams = mfem.MatrixFreeAMS(a, A, fespace_nd, muinv, delta, None, ess_bdr)
        cg = mfem.CGSolver(MPI.COMM_WORLD)
        cg.SetRelTol(1e-12)
        cg.SetMaxIter(1000)
        cg.SetPrintLevel(1)
        cg.SetOperator(*A)
        cg.SetPreconditioner(ams)
        cg.Mult(B, X)

    else:
        AM = A.AsHypreParMatrix()
        if myid == 0:
            print("\nSize of linear system: "+str(AM.GetGlobalNumRows()))

        prec_fespace = a.SCParFESpace() if a.StaticCondensationIsEnabled() else fespace_nd

        ams = mfem.HypreAMS(AM, prec_fespace)
        pcg = mfem.HyprePCG(AM)
        pcg.SetTol(1e-12)
        pcg.SetMaxIter(500)
        pcg.SetPrintLevel(2)
        pcg.SetPreconditioner(ams)
        pcg.Mult(B, X)

    # 15. Recover the parallel grid function corresponding to X. This is the
    #     local finite element solution on each processor.
    a.RecoverFEMSolution(X, b, x)

    # 16. Save the refined mesh and the solution in parallel. This output can
    #     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    pmesh.Print("mesh"+smyid, 8)
    x.Save("sol"+smyid, 8)

    # 17. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock << "solution\n" << pmesh << x
        sol_sock << "window_title 'Vector Potential'"
        sol_sock << "window_geometry 800 0 400 350"
        sol_sock.flush()

    # 18. Compute the magnetic flux as the curl of the solution
    curl = mfem.ParDiscreteLinearOperator(fespace_nd, fespace_rt)
    curl.AddDomainInterpolator(mfem.CurlInterpolator())
    curl.Assemble()
    curl.Finalize()

    dx = mfem.ParGridFunction(fespace_rt)
    curl.Mult(x, dx)

    # 19. Save the curl of the solution in parallel. This output can be viewed
    #     later using GLVis: "glvis -np <np> -m mesh -g dsol".
    dx.Save("dsol"+smyid, 8)

    # 19. Send the curl of the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << pmesh << dx
        sol_sock << "window_title 'Magnetic Flux'"
        sol_sock << "window_geometry 1200 0 400 350"
        sol_sock.flush()


def ComputeCurrentDensityOnSubMesh(order, phi0_attr, phi1_attr, jn_zero_attr,
                                   j_cond):
    # Extract the finite element space and mesh on which j_cond is defined
    fes_cond_rt = j_cond.ParFESpace()
    pmesh_cond = fes_cond_rt.GetParMesh()
    dim = pmesh_cond.Dimension()

    # Define a parallel finite element space on the SubMesh. Here we use the H1
    # finite elements for the electrostatic potential.
    fec_h1 = mfem.H1_FECollection(order, dim)
    fes_cond_h1 = mfem.ParFiniteElementSpace(pmesh_cond, fec_h1)

    # Define the conductivity coefficient and the boundaries associated with the
    # fixed potentials phi0 and phi1 which will drive the current.
    sigmaCoef = mfem.ConstantCoefficient(1.0)
    ess_bdr_phi = mfem.intArray(pmesh_cond.bdr_attributes.Max())
    ess_bdr_j = mfem.intArray(pmesh_cond.bdr_attributes.Max())
    ess_bdr_tdof_phi = mfem.intArray()
    ess_bdr_phi.Assign(0)
    ess_bdr_j.Assign(0)

    for i in range(phi0_attr.Size()):
        ess_bdr_phi[phi0_attr[i]-1] = 1

    for i in range(phi1_attr.Size()):
        ess_bdr_phi[phi1_attr[i]-1] = 1

    for i in range(jn_zero_attr.Size()):
        ess_bdr_j[jn_zero_attr[i]-1] = 1

    fes_cond_h1.GetEssentialTrueDofs(ess_bdr_phi, ess_bdr_tdof_phi)

    # Setup the bilinear form corresponding to -Div(sigma Grad phi)
    a_h1 = mfem.ParBilinearForm(fes_cond_h1)
    a_h1.AddDomainIntegrator(mfem.DiffusionIntegrator(sigmaCoef))
    a_h1.Assemble()

    # Set the r.h.s. to zero
    b_h1 = mfem.ParLinearForm(fes_cond_h1)
    b_h1.Assign(0.0)

    # Setup the boundary conditions on phi
    one = mfem.ConstantCoefficient(1.0)
    zero = mfem.ConstantCoefficient(0.0)
    phi_h1 = mfem.GridFunction(fes_cond_h1)
    phi_h1.Assign(0.0)

    bdr0 = mfem.intArray([0]*pmesh_cond.bdr_attributes.Max())
    for i in range(phi0_attr.Size()):
        bdr0[phi0_attr[i]-1] = 1

    phi_h1.ProjectBdrCoefficient(zero, bdr0)

    bdr1 = mfem.intArray([0]*pmesh_cond.bdr_attributes.Max())
    for i in range(phi1_attr.Size()):
        bdr1[phi1_attr[i]-1] = 1

    phi_h1.ProjectBdrCoefficient(one, bdr1)

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a_h1.FormLinearSystem(ess_bdr_tdof_phi, phi_h1, b_h1, A, X, B)

    # Solve the linear system using algebraic multigrid
    if myid == 0:
        print("\nSolving for electric potential using CG with AMG")
    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a_h1.FormLinearSystem(ess_bdr_tdof_phi, phi_h1, b_h1, A, X, B)

    prec = mfem.HypreBoomerAMG()
    cg = mfem.CGSolver(MPI.COMM_WORLD)
    cg.SetRelTol(1e-12)
    cg.SetMaxIter(2000)
    cg.SetPrintLevel(1)
    cg.SetPreconditioner(prec)
    cg.SetOperator(A.Ptr())
    cg.Mult(B, X)
    a_h1.RecoverFEMSolution(X, b_h1, phi_h1)

    # visualization
    if visualization:
        port_sock = mfem.socketstream("localhost", 19916)
        port_sock << "parallel " << num_procs << " " << myid << "\n"
        port_sock.precision(8)
        port_sock << "solution\n" << pmesh_cond << phi_h1
        port_sock << "window_title 'Conductor Potential'"
        port_sock << "window_geometry 0 0 400 350"
        port_sock.flush()

    # Solve for the current density J = -sigma Grad phi with boundary conditions
    # J.n = 0 on the walls of the conductor but not on the ports where phi=0 and
    # phi=1.

    # J will be computed in H(div) so we need an RT mass matrix
    m_rt = mfem.ParBilinearForm(fes_cond_rt)
    m_rt.AddDomainIntegrator(mfem.VectorFEMassIntegrator())
    m_rt.Assemble()

    # Assemble the (sigma Grad phi) operator
    d_h1 = mfem.ParMixedBilinearForm(fes_cond_h1, fes_cond_rt)
    d_h1.AddDomainIntegrator(mfem.MixedVectorGradientIntegrator(sigmaCoef))
    d_h1.Assemble()

    # Compute the r.h.s, b_rt = sigma E = -sigma Grad phi
    b_rt = mfem.ParLinearForm(fes_cond_rt)
    d_h1.Mult(phi_h1, b_rt)
    b_rt *= -1.0

    # Apply the necessary boundary conditions and solve for J in H(div)
    glb_size_rt = fes_cond_rt.GlobalTrueVSize()
    if myid == 0:
        print("\nSolving for current density in H(Div), using diagonally scaled CG")
        print("Size of linear system: " + str(glb_size_rt))

    ess_bdr_tdof_rt = mfem.intArray()
    M = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()

    fes_cond_rt.GetEssentialTrueDofs(ess_bdr_j, ess_bdr_tdof_rt)
    j_cond.Assign(0.0)
    m_rt.FormLinearSystem(ess_bdr_tdof_rt, j_cond, b_rt, M, X, B)

    pred = mfem.HypreDiagScale()
    cg = mfem.CGSolver(MPI.COMM_WORLD)
    cg.SetRelTol(1e-12)
    cg.SetMaxIter(2000)
    cg.SetPrintLevel(1)
    cg.SetPreconditioner(prec)
    cg.SetOperator(M.Ptr())
    cg.Mult(B, X)
    m_rt.RecoverFEMSolution(X, b_rt, j_cond)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex34p (Sourse Function using a Submehs Transfer )')

    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly in serial.")
    parser.add_argument('-rp', '--refine-parallel',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly in parallel.")
    parser.add_argument("-o", "--order",
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).")
    parser.add_argument("-mc", "--magnetic-cond",
                        action='store', default=1e-6, type=float,
                        help="Magnetic Conductivity")
    parser.add_argument("-sc", "--static-condensation",
                        action='store_true',
                        help="Enable static condensation.")
    parser.add_argument("-hex", "--hex-mesh",
                        action='store_true',  default=False,
                        help="Mixed mesh of hexahedral mesh.")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    parser.add_argument("-no-vis", "--no-visualization",
                        action='store_true', default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    globals()["pa"] = args.partial_assembly
    globals()["visualization"] = not args.no_visualization

    run(ser_ref_levels=args.refine_serial,
        par_ref_levels=args.refine_parallel,
        order=args.order,
        delta_const=args.magnetic_cond,
        static_cond=args.static_condensation,
        mixed=(not args.hex_mesh),
        device_config=args.device)
