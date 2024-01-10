'''
   PyMFEM example 34

   See c++ version in the MFEM library for more detail

   Sample runs:  python ex34.py -o 2
                 python ex34.py -pa -hex

   Device sample runs:
                 python ex34.py -o 2 -pa -hex -d cuda
                 python ex34.py -o 2 -no-pa -d cuda

'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor

pa = False
algebraic_ceed = False
visualization = True


def run(ref_levels=1,
        order=1,
        delta_const=1e-6,
        static_cond=False,
        mixed=True,
        device_config="cpu",
        use_umfpack=False):

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

    # 2. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device_config)
    device.Print()

    # 3. Read the (serial) mesh from the given mesh file on all processors.  We
    #    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    #    and volume meshes with the same code.
    mesh_file = expanduser(
        join(os.path.dirname(__file__), *mesh_file.split("/")))
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()
    if not mixed or pa:
        mesh.UniformRefinement()
        if ref_levels > 0:
            ref_levels = ref_levels - 1
    submesh_attr = -1

    if cond_attr.Size() == 0 and submesh_elems.Size() > 0:
        max_attr = mesh.attributes.Max()
        submesh_attr = max_attr + 1

        for i in range(submesh_elems.Size()):
            mesh.SetAttribute(submesh_elems[i], submesh_attr)
        mesh.SetAttributes()
        if cond_attr.Size() == 0:
            cond_attr.Append(submesh_attr)

    # 4. Refine the serial mesh on all processors to increase the resolution. In
    #    this example we do 'ref_levels' of uniform refinement.
    for l in range(ref_levels):
        mesh.UniformRefinement()

    # 5b. Extract a submesh covering a portion of the domain
    mesh_cond = mfem.SubMesh.CreateFromDomain(mesh, cond_attr)

    # 6. Define a suitable finite element space on the SubMesh and compute
    #    the current density as an H(div) field.
    fec_cond_rt = mfem.RT_FECollection(order - 1, dim)
    fes_cond_rt = mfem.FiniteElementSpace(mesh_cond, fec_cond_rt)
    j_cond = mfem.GridFunction(fes_cond_rt)

    ComputeCurrentDensityOnSubMesh(order, phi0_attr, phi1_attr, jn_zero_attr,
                                   j_cond, use_umfpack)

    # 6a. Save the SubMesh and associated current density in parallel. This
    #     output can be viewed later using GLVis:
    #        "glvis -np <np> -m cond_mesh -g cond_j"
    mesh_cond.Print("cond.mesh", 8)
    j_cond.Save("cond_j.gf", 8)

    # 6b. Send the current density, computed on the SubMesh, to a GLVis server.
    if visualization:
        port_sock = mfem.socketstream("localhost", 19916)
        port_sock.precision(8)
        port_sock << "solution\n" << mesh_cond << j_cond
        port_sock << "window_title 'Conductor J'"
        port_sock << "window_geometry 400 0 400 350"
        port_sock.flush()

    # 7. Define a parallel finite element space on the full mesh. Here we use
    #    the H(curl) finite elements for the vector potential and H(div) for the
    #    current density.
    fec_nd = mfem.ND_FECollection(order, dim)
    fec_rt = mfem.RT_FECollection(order - 1, dim)
    fespace_nd = mfem.FiniteElementSpace(mesh, fec_nd)
    fespace_rt = mfem.FiniteElementSpace(mesh, fec_rt)

    j_full = mfem.GridFunction(fespace_rt)
    j_full.Assign(0.0)
    mesh_cond.Transfer(j_cond, j_full)

    # 7a. Send the transferred current density to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << j_full
        sol_sock << "window_title 'J Full'"
        sol_sock << "window_geometry 400 430 400 350"
        sol_sock.flush()

    # 8. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, the boundary conditions are defined by
    #    marking all the boundary attributes except for those on a symmetry
    #    plane as essential (Dirichlet) and converting them to a list of true
    #    dofs.
    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr.SetSize(mesh.bdr_attributes.Max())
        ess_bdr.Assign(1)
        for i in range(sym_plane_attr.Size()):
            ess_bdr[sym_plane_attr[i]-1] = 0
        fespace_nd.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    #  9. Set up the parallel linear form b(.) which corresponds to the
    #     right-hand side of the FEM linear system, which in this case is
    #     (J,W_i) where J is given by the function H(div) field transferred
    #     from the SubMesh and W_i are the basis functions in the finite
    #     element fespace.
    jCoef = mfem.VectorGridFunctionCoefficient(j_full)
    b = mfem.LinearForm(fespace_nd)
    b.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(jCoef))
    b.Assemble()

    # 10. Define the solution vector x as a parallel finite element grid
    #     function corresponding to fespace. Initialize x to zero.
    x = mfem.GridFunction(fespace_nd)
    x.Assign(0.0)

    # 11. Set up the parallel bilinear form corresponding to the EM diffusion
    #     operator curl muinv curl + delta I, by adding the curl-curl and the
    #     mass domain integrators. For standard magnetostatics equations choose
    #     delta << 1. Larger values of delta should make the linear system
    #     easier to solve at the expense of resembling a diffusive quasistatic
    #     magnetic field.  A reasonable balance must be found whenever the mesh
    #     or problem setup is altered.
    muinv = mfem.ConstantCoefficient(1.0)
    delta = mfem.ConstantCoefficient(delta_const)
    a = mfem.BilinearForm(fespace_nd)
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(delta))

    # 12. Assemble the parallel bilinear form and the corresponding linear
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

    # 13. Solve the system AX=B
    if pa:  # Jacobi preconditioning in partial assembly mode
        print("\n".join(("\nSolving for magnetic vector potential ",
                         "using CG with a Jacobi preconditioner")))

        M = mfem.OperatorJacobiSmoother(a, ess_tdof_list)
        mfem.PCG(A, M, B, X, 1, 1000, 1e-12, 0.0)

    else:
        if not use_umfpack:
            print("\n".join(("\nSolving for magnetic vector potential ",
                             "using CG with a Gauss-Seidel preconditioner")))

            # 13a. Define a simple symmetric Gauss-Seidel preconditioner and use
            #     it to solve the system Ax=b with PCG.
            AA = mfem.OperatorHandle2SparseMatrix(A)
            M = mfem.GSSmoother(AA)
            mfem.PCG(A, M, B, X, 1, 500, 1e-12, 0.0)

        else:
            print("".join(("\nSolving for magnetic vector potential ",
                           "using UMFPack")))

            # 13a. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the
            #     system.
            umf_solver = mfem.UMFPackSolver()
            umf_solver.Control[
                mfem.UMFPACK_ORDERING] = mfem.UMFPACK_ORDERING_METIS
            umf_solver.SetOperator(A)
            umf_solver.Mult(B, X)

    # 14. Recover the parallel grid function corresponding to X. This is the
    #     local finite element solution on each processor.
    a.RecoverFEMSolution(X, b, x)

    # 15. Save the refined mesh and the solution in parallel. This output can
    #     be viewed later using GLVis: "glvis -m refined.mesh -g sol.gf".
    mesh.Print("refined.mesh", 8)
    x.Save("sol.gf", 8)

    # 16. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "solution\n" << mesh << x
        sol_sock << "window_title 'Vector Potential'"
        sol_sock << "window_geometry 800 0 400 350"
        sol_sock.flush()

    # 17. Compute the magnetic flux as the curl of the solution
    curl = mfem.DiscreteLinearOperator(fespace_nd, fespace_rt)
    curl.AddDomainInterpolator(mfem.CurlInterpolator())
    curl.Assemble()
    curl.Finalize()

    dx = mfem.GridFunction(fespace_rt)
    curl.Mult(x, dx)

    # 18. Save the curl of the solution in parallel. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g dsol.gf".
    dx.Save("dsol.gf", 8)

    # 19. Send the curl of the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << dx
        sol_sock << "window_title 'Magnetic Flux'"
        sol_sock << "window_geometry 1200 0 400 350"
        sol_sock.flush()


def ComputeCurrentDensityOnSubMesh(order, phi0_attr, phi1_attr, jn_zero_attr,
                                   j_cond, use_umfpack):
    # Extract the finite element space and mesh on which j_cond is defined
    fes_cond_rt = j_cond.FESpace()
    mesh_cond = fes_cond_rt.GetMesh()
    dim = mesh_cond.Dimension()

    # Define a parallel finite element space on the SubMesh. Here we use the H1
    # finite elements for the electrostatic potential.
    fec_h1 = mfem.H1_FECollection(order, dim)
    fes_cond_h1 = mfem.FiniteElementSpace(mesh_cond, fec_h1)

    # Define the conductivity coefficient and the boundaries associated with the
    # fixed potentials phi0 and phi1 which will drive the current.
    sigmaCoef = mfem.ConstantCoefficient(1.0)
    ess_bdr_phi = mfem.intArray(mesh_cond.bdr_attributes.Max())
    ess_bdr_j = mfem.intArray(mesh_cond.bdr_attributes.Max())
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
    a_h1 = mfem.BilinearForm(fes_cond_h1)
    a_h1.AddDomainIntegrator(mfem.DiffusionIntegrator(sigmaCoef))
    a_h1.Assemble()

    # Set the r.h.s. to zero
    b_h1 = mfem.LinearForm(fes_cond_h1)
    b_h1.Assign(0.0)

    # Setup the boundary conditions on phi
    one = mfem.ConstantCoefficient(1.0)
    zero = mfem.ConstantCoefficient(0.0)
    phi_h1 = mfem.GridFunction(fes_cond_h1)
    phi_h1.Assign(0.0)

    bdr0 = mfem.intArray([0]*mesh_cond.bdr_attributes.Max())
    for i in range(phi0_attr.Size()):
        bdr0[phi0_attr[i]-1] = 1

    phi_h1.ProjectBdrCoefficient(zero, bdr0)

    bdr1 = mfem.intArray([0]*mesh_cond.bdr_attributes.Max())
    for i in range(phi1_attr.Size()):
        bdr1[phi1_attr[i]-1] = 1

    phi_h1.ProjectBdrCoefficient(one, bdr1)

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a_h1.FormLinearSystem(ess_bdr_tdof_phi, phi_h1, b_h1, A, X, B)

    #  Solve the linear system

    if not pa:
        if not use_umfpack:
            print("".join(("\nSolving for electric potential using PCG ",
                           "with a Gauss-Seidel preconditioner")))

            # Use a simple symmetric Gauss-Seidel preconditioner with PCG.
            AA = mfem.OperatorHandle2SparseMatrix(A)
            M = mfem.GSSmoother(AA)
            mfem.PCG(AA, M, B, X, 1, 200, 1e-12, 0.0)
        else:
            print("\nSolving for electric potential using UMFPack")

            #  If MFEM was compiled with SuiteSparse,
            # use UMFPACK to solve the system.
            umf_solver = mfem.UMFPackSolver()
            umf_solver.Control[mfem.UMFPACK_ORDERING] = mfem.UMFPACK_ORDERING_METIS
            umf_solver.SetOperator(A)
            umf_solver.Mult(B, X)
    else:
        print("\nSolving for electric potential using CG")

        if mfem.UsesTensorBasis(fes_cond_h1):
            if algebraic_ceed:
                assert False, "not supported"
                # prec = new ceed::AlgebraicSolver(a, ess_tdof_list);
                # PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
            else:
                M = mfem.OperatorJacobiSmoother(a_h1, ess_bdr_tdof_phi)
                mfem.PCG(A, M, B, X, 1, 400, 1e-12, 0.0)
        else:
            mfem.CG(A, B, X, 1, 400, 1e-12, 0.0)
    a_h1.RecoverFEMSolution(X, b_h1, phi_h1)

    # visualization
    if visualization:
        port_sock = mfem.socketstream("localhost", 19916)
        port_sock.precision(8)
        port_sock << "solution\n" << mesh_cond << phi_h1
        port_sock << "window_title 'Conductor Potential'"
        port_sock << "window_geometry 0 0 400 350"
        port_sock.flush()

    # Solve for the current density J = -sigma Grad phi with boundary conditions
    # J.n = 0 on the walls of the conductor but not on the ports where phi=0 and
    # phi=1.

    # J will be computed in H(div) so we need an RT mass matrix
    m_rt = mfem.BilinearForm(fes_cond_rt)
    m_rt.AddDomainIntegrator(mfem.VectorFEMassIntegrator())
    m_rt.Assemble()

    # Assemble the (sigma Grad phi) operator
    d_h1 = mfem.MixedBilinearForm(fes_cond_h1, fes_cond_rt)
    d_h1.AddDomainIntegrator(mfem.MixedVectorGradientIntegrator(sigmaCoef))
    d_h1.Assemble()

    # Compute the r.h.s, b_rt = sigma E = -sigma Grad phi
    b_rt = mfem.LinearForm(fes_cond_rt)
    d_h1.Mult(phi_h1, b_rt)
    b_rt *= -1.0

    # Apply the necessary boundary conditions and solve for J in H(div)
    print("\nSolving for current density in H(Div), using diagonally scaled CG")
    print("Size of linear system: " + str(fes_cond_rt.GetTrueVSize()))

    ess_bdr_tdof_rt = mfem.intArray()
    M = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()

    fes_cond_rt.GetEssentialTrueDofs(ess_bdr_j, ess_bdr_tdof_rt)
    j_cond.Assign(0.0)
    m_rt.FormLinearSystem(ess_bdr_tdof_rt, j_cond, b_rt, M, X, B)

    cg = mfem.CGSolver()
    cg.SetRelTol(1e-12)
    cg.SetMaxIter(2000)
    cg.SetPrintLevel(1)
    cg.SetOperator(M)
    cg.Mult(B, X)
    m_rt.RecoverFEMSolution(X, b_rt, j_cond)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex34 (Sourse Function using a Submehs Transfer )')

    parser.add_argument('-r', '--refine',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly.")
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
                        action='store_false', default=True,
                        help='Enable GLVis visualization')

    if hasattr(mfem, "UMFPackSolver"):
        parser.add_argument("-noumfpack", "--no-use-umfpack",
                            action='store_true', default=False,
                            help='Enable UMFPACK')

    args = parser.parse_args()
    parser.print_options(args)

    globals()["pa"] = args.partial_assembly
    globals()["visualization"] = args.no_visualization

    if hasattr(mfem, "UMFPackSolver"):
        use_umfpack = True if not args.no_use_umfpack else False
    else:
        use_umfpack = False

    run(ref_levels=args.refine,
        order=args.order,
        delta_const=args.magnetic_cond,
        static_cond=args.static_condensation,
        mixed=(not args.hex_mesh),
        device_config=args.device,
        use_umfpack=use_umfpack)
