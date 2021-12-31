'''
   MFEM example 22
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi


def run(mesh_file="",
        order=1,
        ref_levels=0,
        visualization=1,
        prob=0,
        harm_conv=True,
        a_coef=0.0,
        epsilon=1.0,
        sigma=20.,
        mu=1.0,
        device='cpu',
        pa=False,
        freq=-1.0,
        numba=False):

    if a_coef != 0.0:
        mu = 1.0 / a_coef
    omega = 10.
    if freq > 0.0:
        omega = 2.0 * pi * freq
    exact_sol = check_for_inline_mesh(mesh_file)
    if exact_sol:
        print("Identified a mesh with known exact solution")

    conv = mfem.ComplexOperator.HERMITIAN if harm_conv else mfem.ComplexOperator.BLOCK_SYMMETRIC

    # 2. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device)
    device.Print()

    # 3. Read the mesh from the given mesh file. We can handle triangular,
    #    quadrilateral, tetrahedral, hexahedral, surface and volume meshes
    #    with the same code.
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()

    # 4. Refine the mesh to increase resolution. In this example we do
    #    'ref_levels' of uniform refinement where the user specifies
    #    the number of levels with the '-r' option.
    for i in range(ref_levels):
        mesh.UniformRefinement()

    # 5. Define a finite element space on the mesh. Here we use continuous
    #    Lagrange, Nedelec, or Raviart-Thomas finite elements of the specified
    #    order.
    if (dim == 1 and prob != 0):
        print("Switching to problem type 0, H1 basis functions, for 1 dimensional mesh.")
        prob = 0

    if prob == 0:
        fec = mfem.H1_FECollection(order, dim)
    elif prob == 1:
        fec = mfem.ND_FECollection(order, dim)
    elif prob == 2:
        fec = mfem.RT_FECollection(order-1, dim)
    else:
        assert False, "unknown problem"

    fespace = mfem.FiniteElementSpace(mesh, fec)
    print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))

    # 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    #    In this example, the boundary conditions are defined based on the type
    #    of mesh and the problem type.
    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray()

    if mesh.bdr_attributes.Size():
        ess_bdr.SetSize(mesh.bdr_attributes.Max())
        ess_bdr.Assign(1)
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 7. Set up the linear form b(.) which corresponds to the right-hand side of
    #    the FEM linear system.
    b = mfem.ComplexLinearForm(fespace, conv)
    b.Assign(0.0)

    # 8. Define the solution vector u as a complex finite element grid function
    #    corresponding to fespace. Initialize u with initial guess of 1+0i or
    #    the exact solution if it is known.

    # 8-a. Prepare coefficient in pure Python and numba
    sdim = mesh.SpaceDimension()

    if numba:
        @mfem.jit.scalar()
        def u0_real_exact(x):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            return exp(-1j * kappa * x[dim - 1]).real

        @mfem.jit.scalar()
        def u0_imag_exact(x):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            return exp(-1j * kappa * x[dim - 1]).imag

        @mfem.jit.vector()
        def u1_real_exact(x, out):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            out[0] = exp(-1j * kappa * x[dim - 1]).real
            out[1] = 0
            out[2] = 0

        @mfem.jit.vector()
        def u1_imag_exact(x, out):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            out[0] = exp(-1j * kappa * x[dim - 1]).imag
            out[1] = 0
            out[2] = 0

        @mfem.jit.vector()
        def u2_real_exact(x, out):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            for i in range(dim):
                out[i] = 0
            out[dim-1] = exp(-1j * kappa * x[dim - 1]).real

        @mfem.jit.vector()
        def u2_imag_exact(x, out):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            for i in range(dim):
                out[i] = 0
            out[dim-1] = exp(-1j * kappa * x[dim - 1]).imag

    else:
        def u0_exact(x):
            alpha = (epsilon * omega - 1j * sigma)
            kappa = sqrt(mu * omega * alpha)
            return exp(-1j * kappa * x[dim - 1])

        class cu0_real_exact(mfem.PyCoefficient):
            def EvalValue(self, x):
                return u0_exact(x).real

        class cu0_imag_exact(mfem.PyCoefficient):
            def EvalValue(self, x):
                return u0_exact(x).imag
        u0_real_exact = cu0_real_exact()
        u0_imag_exact = cu0_imag_exact()

        class cu1_real_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                ret = x*0
                ret[0] = u0_exact(x).real
                return ret

        class cu1_imag_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                ret = x*0
                ret[0] = u0_exact(x).imag
                return ret

        u1_real_exact = cu1_real_exact()
        u1_imag_exact = cu1_imag_exact()

        class cu2_real_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                ret = x*0
                ret[dim-1] = u0_exact(x).real
                return ret

        class cu2_imag_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                ret = x*0
                ret[dim-1] = u0_exact(x).imag
                return ret

        u2_real_exact = cu2_real_exact()
        u2_imag_exact = cu2_imag_exact()

    u = mfem.ComplexGridFunction(fespace)
    if exact_sol:
        u_exact = mfem.ComplexGridFunction(fespace)

    u0_r = u0_real_exact
    u0_i = u0_imag_exact
    u1_r = u1_real_exact
    u1_i = u1_imag_exact
    u2_r = u2_real_exact
    u2_i = u2_imag_exact

    zeroCoef = mfem.ConstantCoefficient(0.0)
    oneCoef = mfem.ConstantCoefficient(1.0)

    zeroVec = mfem.Vector(dim)
    zeroVec.Assign(0.0)
    oneVec = mfem.Vector(dim)
    oneVec.Assign(0.0)
    oneVec[dim-1 if prob == 2 else 0] = 1.0
    zeroVecCoeff = mfem.VectorConstantCoefficient(zeroVec)
    oneVecCoeff = mfem.VectorConstantCoefficient(oneVec)

    if prob == 0:
        if exact_sol:
            u.ProjectBdrCoefficient(u0_r, u0_i, ess_bdr)
            u_exact.ProjectCoefficient(u0_r, u0_i)
        else:
            u.ProjectBdrCoefficient(oneCoef, zeroCoef, ess_bdr)
    elif prob == 1:
        if exact_sol:
            u.ProjectBdrCoefficientTangent(u1_r, u1_i, ess_bdr)
            u_exact.ProjectCoefficient(u1_r, u1_i)
        else:
            u.ProjectBdrCoefficientTangent(oneVecCoef, zeroVecCoef, ess_bdr)
    elif prob == 2:
        if exact_sol:
            u.ProjectBdrCoefficientNormal(u2_r, u2_i, ess_bdr)
            u_exact.ProjectCoefficient(u2_r, u2_i)
        else:
            u.ProjectBdrCoefficientNormal(oneVecCoef, zeroVecCoef, ess_bdr)

    if (visualization and exact_sol):
        sol_sock_r = mfem.socketstream("localhost", 19916)
        sol_sock_r.precision(8)
        sol_sock_r << "solution\n" << mesh << u_exact.real(
        ) << "window_title 'Exact: Real Part'"
        sol_sock_r.flush()
        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << mesh << u_exact.imag(
        ) << "window_title 'Exact: Imag Part'"
        sol_sock_i.flush()

    # 9. Set up the sesquilinear form a(.,.) on the finite element space
    #    corresponding to the damped harmonic oscillator operator of the
    #    appropriate type:
    #
    #    0) A scalar H1 field
    #       -Div(a Grad) - omega^2 b + i omega c
    #
    #    1) A vector H(Curl) field
    #       Curl(a Curl) - omega^2 b + i omega c
    #
    #    2) A vector H(Div) field
    #       -Grad(a Div) - omega^2 b + i omega c

    stiffnessCoef = mfem.ConstantCoefficient(1.0/mu)
    massCoef = mfem.ConstantCoefficient(-omega * omega * epsilon)
    lossCoef = mfem.ConstantCoefficient(omega * sigma)
    negMassCoef = mfem.ConstantCoefficient(omega * omega * epsilon)

    a = mfem.SesquilinearForm(fespace, conv)
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    if prob == 0:
        a.AddDomainIntegrator(mfem.DiffusionIntegrator(stiffnessCoef),
                              None)
        a.AddDomainIntegrator(mfem.MassIntegrator(massCoef),
                              mfem.MassIntegrator(lossCoef))
    elif prob == 1:
        a.AddDomainIntegrator(mfem.CurlCurlIntegrator(stiffnessCoef),
                              None)
        a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(massCoef),
                              mfem.VectorFEMassIntegrator(lossCoef))
    elif prob == 2:
        a.AddDomainIntegrator(mfem.DivDivIntegrator(stiffnessCoef),
                              None)
        a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(massCoef),
                              mfem.VectorFEMassIntegrator(lossCoef))
    else:
        assert False, "Unknown probm"

    # 9a. Set up the bilinear form for the preconditioner corresponding to the
    #     appropriate operator
    #
    #      0) A scalar H1 field
    #         -Div(a Grad) - omega^2 b + omega c
    #
    #      1) A vector H(Curl) field
    #         Curl(a Curl) + omega^2 b + omega c
    #
    #      2) A vector H(Div) field
    #         -Grad(a Div) - omega^2 b + omega c

    pcOp = mfem.BilinearForm(fespace)
    if pa:
        pcOp.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

    if prob == 0:
        pcOp.AddDomainIntegrator(mfem.DiffusionIntegrator(stiffnessCoef))
        pcOp.AddDomainIntegrator(mfem.MassIntegrator(massCoef))
        pcOp.AddDomainIntegrator(mfem.MassIntegrator(lossCoef))
    elif prob == 1:
        pcOp.AddDomainIntegrator(mfem.CurlCurlIntegrator(stiffnessCoef))
        pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(negMassCoef))
        pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(lossCoef))
    elif prob == 2:
        pcOp.AddDomainIntegrator(mfem.DivDivIntegrator(stiffnessCoef))
        pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(massCoef))
        pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(lossCoef))

    # 10. Assemble the form and the corresponding linear system, applying any
    #     necessary transformations such as: assembly, eliminating boundary
    #     conditions, conforming constraints for non-conforming AMR, etc.
    a.Assemble()
    pcOp.Assemble()

    A = mfem.OperatorHandle()
    B = mfem.Vector()
    U = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, u, b, A, U, B)

    print("Size of linear system: " + str(A.Width()))

    # 11. Define and apply a GMRES solver for AU=B with a block diagonal
    #     preconditioner based on the appropriate sparse smoother.
    blockOffsets = mfem.intArray()
    blockOffsets.SetSize(3)
    blockOffsets[0] = 0
    blockOffsets[1] = A.Height() // 2
    blockOffsets[2] = A.Height() // 2
    blockOffsets.PartialSum()

    BDP = mfem.BlockDiagonalPreconditioner(blockOffsets)

    if pa:
        pc_r = mfem.OperatorJacobiSmoother(pcOp, ess_tdof_list)
        pc_i = None
    else:
        PCOp = mfem.OperatorHandle()
        pcOp.SetDiagonalPolicy(mfem.Operator.DIAG_ONE)
        pcOp.FormSystemMatrix(ess_tdof_list, PCOp)
        if prob == 0:
            pc_r = mfem.DSmoother(mfem.OperatorHandle2SparseMatrix(PCOp))
        elif prob == 1:
            pc_r = mfem.GSSmoother(mfem.OperatorHandle2SparseMatrix(PCOp))
        elif prob == 2:
            pc_r = mfem.DSmoother(mfem.OperatorHandle2SparseMatrix(PCOp))

    s = 1.0 if prob != 1 else -1.0

    pc_i = mfem.ScaledOperator(pc_r,
                               s if conv == mfem.ComplexOperator.HERMITIAN else -s)
    BDP.SetDiagonalBlock(0, pc_r)
    BDP.SetDiagonalBlock(1, pc_i)
    # BDP.owns_blocks = 1  #this is not necessary since Python owns smoothers

    gmres = mfem.GMRESSolver()
    gmres.SetPreconditioner(BDP)
    gmres.SetOperator(A.Ptr())
    gmres.SetRelTol(1e-12)
    gmres.SetMaxIter(1000)
    gmres.SetPrintLevel(1)
    gmres.Mult(B, U)

    # 12. Recover the solution as a finite element grid function and compute the
    #     errors if the exact solution is known.
    a.RecoverFEMSolution(U, b, u)

    if exact_sol:
        err_r = -1.0
        err_i = -1.0

        if prob == 0:
            err_r = u.real().ComputeL2Error(u0_r)
            err_i = u.imag().ComputeL2Error(u0_i)
        elif prob == 1:
            err_r = u.real().ComputeL2Error(u1_r)
            err_i = u.imag().ComputeL2Error(u1_i)
        elif prob == 2:
            err_r = u.real().ComputeL2Error(u2_r)
            err_i = u.imag().ComputeL2Error(u2_i)

    print()
    print("|| Re (u_h - u) ||_{L^2} = " + "{:g}".format(err_r))
    print("|| Im (u_h - u) ||_{L^2} = " + "{:g}".format(err_i))

    # 13. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m mesh -g sol".
    mesh.Print("refined.mesh", 8)
    u.real().Save("sol_r.gf", 8)
    u.imag().Save("sol_i.gf", 8)

    # 14. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock_r = mfem.socketstream("localhost", 19916)
        sol_sock_r.precision(8)
        sol_sock_r << "solution\n" << mesh << u.real(
        ) << "window_title 'Solution: Real Part'"
        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << mesh << u.imag(
        ) << "window_title 'Solution: Imag Part'"

    if visualization and exact_sol:
        u_exact -= u
        sol_sock_r = mfem.socketstream("localhost", 19916)
        sol_sock_r.precision(8)
        sol_sock_r << "solution\n" << mesh << u_exact.real(
        ) << "window_title 'Error: Real Part'"
        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << mesh << u_exact.imag(
        ) << "window_title 'Error: Imag Part'"

    if visualization:
        u_t = mfem.GridFunction(fespace)
        u_t.Assign(u.real())
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << u_t
        sol_sock << "window_title 'Harmonic Solution (t = 0.0 T)'"
        sol_sock << "pause\n"
        sol_sock.flush()

        print("GLVis visualization paused. Press space (in the GLVis window) to resume it.")
        num_frames = 32
        i = 0

        # Let's plot one wave cycle...
        for i in range(32):
            t = (i % num_frames) / num_frames
            oss = "Harmonic Solution (t = " + str(t) + " T)"
            dd = (cos(2.0 * pi * t)*u.real().GetDataArray() +
                  sin(2.0 * pi * t)*u.imag().GetDataArray())
            u_t.Assign(mfem.Vector())
            sol_sock << "solution\n" << mesh << u_t
            sol_sock << "window_title '" << oss << "'"
            sol_sock.flush()


def check_for_inline_mesh(mesh_file):
    return os.path.basename(mesh_file)[:7] == "inline-"


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex22 (Complex linear systems)')
    parser.add_argument('-m', '--mesh',
                        default='inline-quad.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-p",
                        "--problem-type",
                        action='store', type=int, default=0,
                        help="\n".join(["Choose between 0: H_1, 1: H(Curl), or 2: H(Div) "
                                        "damped harmonic oscillator."]))
    parser.add_argument("-r", "--refine",
                        action='store', type=int, default=0,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument("-a", "--stiffness-coef",
                        action='store', type=float, default=0.0,
                        help="Stiffness coefficient (spring constant or 1/mu).")
    parser.add_argument("-mu", "--permeability",
                        action='store', type=float, default=1.0,
                        help="Permeability of free space (or 1/(spring constant)).")
    parser.add_argument("-eps", "--permittivity",
                        action='store', type=float, default=1.0,
                        help="Permittivity of free space (or mass constant).")
    parser.add_argument("-sigma", "--conductivity",
                        action='store', type=float, default=20.0,
                        help="Conductivity (or damping constant).")
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=-1.0,
                        help="Set the frequency for the exact")
    parser.add_argument("-no-herm", "--no-hermitian",
                        action='store_false',
                        default=True,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=True,
                        help='Enable GLVis visualization')
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")

    try:
        from numba import jit
        HAS_NUMBA = True
    except ImportError:
        HAS_NUMBA = False
    parser.add_argument("-n", "--numba",
                        default=int(HAS_NUMBA),
                        type=int,
                        action='store',
                        help="Use Number compiled coefficient")

    args = parser.parse_args()
    args.numba = bool(args.numba)
    parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    run(mesh_file=meshfile,
        order=args.order,
        ref_levels=args.refine,
        prob=args.problem_type,
        visualization=args.visualization,
        # visualization=False,
        harm_conv=args.no_hermitian,
        a_coef=args.stiffness_coef,
        epsilon=args.permittivity,
        sigma=args.conductivity,
        mu=args.permeability,
        device=args.device,
        pa=args.partial_assembly,
        freq=args.frequency,
        numba=args.numba)
