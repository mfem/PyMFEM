'''
   MFEM example 22p
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(mesh_file="",
        order=1,
        ser_ref_levels=0,
        par_ref_levels=0,
        visualization=1,
        prob=0,
        herm_conv=True,
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

    conv = mfem.ComplexOperator.HERMITIAN if herm_conv else mfem.ComplexOperator.BLOCK_SYMMETRIC

    # 2. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device)
    if myid == 0:
        device.Print()    

    # 3. Read the mesh from the given mesh file. We can handle triangular,
    #    quadrilateral, tetrahedral, hexahedral, surface and volume meshes
    #    with the same code.
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()

    # 4. Refine the mesh to increase resolution. In this example we do
    #    'ref_levels' of uniform refinement where the user specifies
    #    the number of levels with the '-r' option.
    for i in range(ser_ref_levels):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    for i in range(par_ref_levels):
        pmesh.UniformRefinement()

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

    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    global_dof = fespace.GlobalTrueVSize()
    if myid == 0:
        print("Number of finite element unknowns: " + str(global_dof))

    # 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    #    In this example, the boundary conditions are defined based on the type
    #    of mesh and the problem type.
    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray()

    if pmesh.bdr_attributes.Size():
        ess_bdr.SetSize(pmesh.bdr_attributes.Max())
        ess_bdr.Assign(1)
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 7. Set up the linear form b(.) which corresponds to the right-hand side of
    #    the FEM linear system.
    b = mfem.ParComplexLinearForm(fespace, conv)
    b.Assign(0.0)

    # 8. Define the solution vector u as a complex finite element grid function
    #    corresponding to fespace. Initialize u with initial guess of 1+0i or
    #    the exact solution if it is known.

    # 8-a. Prepare coefficient in pure Python and numba
    if numba:
        '''
        @mfem.jit.vector()
        def E_exact(x, out):
            out[0] = sin(kappa*x[1])
            out[1] = sin(kappa*x[2])
            out[2] = sin(kappa*x[0])
        @mfem.jit.vector()
        def f_exact(x, out):
            out[0] = (1 + kappa**2)*sin(kappa * x[1])
            out[1] = (1 + kappa**2)*sin(kappa * x[2])
            out[2] = (1 + kappa**2)*sin(kappa * x[0])
        '''
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

        sdim = pmesh.SpaceDimension()

        class cu1_real_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                return (u0_exact(x).real, 0, 0)

        class cu1_imag_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                return (u0_exact(x).imag, 0, 0)
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

    u = mfem.ParComplexGridFunction(fespace)
    if exact_sol:
        u_exact = mfem.ParComplexGridFunction(fespace)

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
        sol_sock_r << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_r << "solution\n" << pmesh << u_exact.real(
        ) << "window_title 'Exact: Real Part'"
        sol_sock_r.flush()
        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << pmesh << u_exact.imag(
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

    a = mfem.ParSesquilinearForm(fespace, conv)
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

    pcOp = mfem.ParBilinearForm(fespace)
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

    if myid == 0:
        print("Size of linear system: " + str(2*global_dof))

    # 11. Define and apply a GMRES solver for AU=B with a block diagonal
    #     preconditioner based on the appropriate sparse smoother.
    blockTrueOffsets = mfem.intArray()
    blockTrueOffsets.SetSize(3)
    blockTrueOffsets[0] = 0
    blockTrueOffsets[1] = A.Height() // 2
    blockTrueOffsets[2] = A.Height() // 2
    blockTrueOffsets.PartialSum()

    BDP = mfem.BlockDiagonalPreconditioner(blockTrueOffsets)

    if pa:
        pc_r = mfem.OperatorJacobiSmoother(pcOp, ess_tdof_list)
        pc_i = None
    else:
        PCOp = mfem.OperatorHandle()
        pcOp.SetDiagonalPolicy(mfem.Operator.DIAG_ONE)
        pcOp.FormSystemMatrix(ess_tdof_list, PCOp)
        if prob == 0:
            pc_r = mfem.HypreBoomerAMG(PCOp.AsHypreParMatrix())
        elif prob == 1:
            pc_r = mfem.HypreAMS(PCOp.AsHypreParMatrix(), fespace)
        elif prob == 2:
            if dim == 2:
                pc_r = mfem.HypreAMS(PCOp.AsHypreParMatrix(), fespace)
            else:
                pc_r = mfem.HypreADS(PCOp.AsHypreParMatrix(), fespace)

    pc_i = mfem.ScaledOperator(pc_r,
                               -1 if conv == mfem.ComplexOperator.HERMITIAN else 1)
    BDP.SetDiagonalBlock(0, pc_r)
    BDP.SetDiagonalBlock(1, pc_i)
    # BDP.owns_blocks = 1  #this is not necessary since Python owns smoothers

    fgmres = mfem.FGMRESSolver(MPI.COMM_WORLD)
    fgmres.SetPreconditioner(BDP)
    fgmres.SetOperator(A.Ptr())
    fgmres.SetRelTol(1e-12)
    fgmres.SetMaxIter(1000)
    fgmres.SetPrintLevel(1)
    fgmres.Mult(B, U)

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

    if myid == 0:
        print()
        print("|| Re (u_h - u) ||_{L^2} = " + "{:g}".format(err_r))
        print("|| Im (u_h - u) ||_{L^2} = " + "{:g}".format(err_i))

    # 13. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m mesh -g sol".
    pmesh.Print("mesh"+smyid, 8)
    u.real().Save("sol_r"+smyid, 8)
    u.imag().Save("sol_i"+smyid, 8)

    # 14. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock_r = mfem.socketstream("localhost", 19916)
        sol_sock_r << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_r.precision(8)
        sol_sock_r << "solution\n" << pmesh << u.real(
        ) << "window_title 'Solution: Real Part'"

        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << pmesh << u.imag(
        ) << "window_title 'Solution: Imag Part'"

    if visualization and exact_sol:
        u_exact -= u
        sol_sock_r = mfem.socketstream("localhost", 19916)
        sol_sock_r << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_r.precision(8)
        sol_sock_r << "solution\n" << pmesh << u_exact.real(
        ) << "window_title 'Error: Real Part'"
        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << pmesh << u_exact.imag(
        ) << "window_title 'Error: Imag Part'"

    if visualization:
        u_t = mfem.ParGridFunction(fespace)
        u_t.Assign(u.real())
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << pmesh << u_t
        sol_sock << "window_title 'Harmonic Solution (t = 0.0 T)'"
        sol_sock << "pause\n"
        sol_sock.flush()

        if myid == 0:
            print(
                "GLVis visualization paused. Press space (in the GLVis window) to resume it.")
        num_frames = 32
        i = 0

        # Let's plot one wave cycle...
        for i in range(32):
            t = (i % num_frames) / num_frames
            oss = "Harmonic Solution (t = " + str(t) + " T)"
            dd = (cos(2.0 * pi * t)*u.real().GetDataArray() +
                  sin(2.0 * pi * t)*u.imag().GetDataArray())
            # we can not load numpy directly...(sorry)
            u_t.Assign(mfem.Vector(dd))
            sol_sock << "parallel " << num_procs << " " << myid << "\n"
            sol_sock << "solution\n" << pmesh << u_t
            sol_sock << "window_title '" << oss << "'"
            sol_sock.flush()


def check_for_inline_mesh(mesh_file):
    return os.path.basename(mesh_file)[:7] == "inline-"


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex22p (Complex linear systems)')
    parser.add_argument('-m', '--mesh',
                        default='inline-quad.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-rs", "--refine-serial",
                        action='store', type=int, default=1,
                        help="Number of times to refine the mesh uniformly in serial.")
    parser.add_argument("-rp", "--refine-parallel",
                        action='store', type=int, default=1,
                        help="Number of times to refine the mesh uniformly in paralle.")
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-p",
                        "--problem-type",
                        action='store', type=int, default=0,
                        help="\n".join(["Choose between 0: H_1, 1: H(Curl), or 2: H(Div) "
                                        "damped harmonic oscillator."]))
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
    parser.add_argument("-herm", "--hermitian",
                        action='store_true',
                        default=True,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument("-no-herm", "--no-hermitian",
                        action='store_true',
                        default=False,
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
    parser.add_argument("-n", "--numba",
                        default=False, type=bool,
                        help="Use Number compiled coefficient")

    args = parser.parse_args()

    herm = False if args.no_hermitian else True
    args.hermitian = herm
    args.no_hermitian = not herm
    if myid == 0:
        parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    run(mesh_file=meshfile,
        order=args.order,
        ser_ref_levels=args.refine_serial,
        par_ref_levels=args.refine_parallel,
        prob=args.problem_type,
        visualization=args.visualization,
        # visualization=False,
        herm_conv=herm,
        a_coef=args.stiffness_coef,
        epsilon=args.permittivity,
        sigma=args.conductivity,
        mu=args.permeability,
        device=args.device,
        pa=args.partial_assembly,
        freq=args.frequency,
        numba=args.numba)
