'''
   MFEM example 26
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(order=1,
        href=1,
        mesh_file='',
        nd=2,
        initref=1,
        length=1.0,
        k=0.5,
        sol=1,
        visualization=True):

    if nd == 2:
        mesh = mfem.Mesh_MakeCartesian2D(
            1, 1, mfem.Element.QUADRILATERAL, True, length, length, False)
    else:
        mesh = mfem.Mesh_MakeCartesian3D(
            1, 1, 1, mfem.Element.HEXAHEDRON, True, length, length, length)

    sdim = mesh.SpaceDimension()
    dim = mesh.Dimension()

    for i in range(initref):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    mesh.Clear()
    pmesh.UniformRefinement()

    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)

    #fespaces = mfem.vector_ParFiniteElementSpace(href+1)
    fespaces = [None]*(href + 1)
    ParMeshes = [None]*(href + 1)
    P = [None] * href  # mfem.vector_HypreParMatrix(href)

    for i in range(href):
        ParMeshes[i] = mfem.ParMesh(pmesh)
        fespaces[i] = mfem.ParFiniteElementSpace(fespace, ParMeshes[i])

        pmesh.UniformRefinement()
        fespace.Update()
        Tr = mfem.OperatorHandle(mfem.Operator.Hypre_ParCSR)
        fespace.GetTrueTransferOperator(fespaces[i], Tr)
        Tr.SetOperatorOwner(False)
        P[i] = Tr.AsHypreParMatrix()

    fespaces[href] = mfem.ParFiniteElementSpace(fespace)

    omega = 2.0 * pi * k

    #
    #  functions (these are JITed using Numba insdie run())
    #

    # local parameters used in JIT function
    from numba import types

    params = {"dim": dim, "omega": omega, "sol": sol}

    sig = types.void(types.CPointer(types.double),
                     types.float64[:],
                     types.float64[:])

    get_sol_Rd = mfem.jit.func(sig, params=params)(get_helmholtz_solution_Re)
    get_sol_Im = mfem.jit.func(sig, params=params)(get_helmholtz_solution_Im)

    params["get_helmholtz_solution_Re"] = get_sol_Rd
    params["get_helmholtz_solution_Im"] = get_sol_Im

    f_Re = mfem.jit.scalar(params=params)(f_exact_Re)
    g_Re = mfem.jit.scalar(params=params)(g_exact_Re)
    grad_Re = mfem.jit.vector(sdim=sdim, params=params)(grad_exact_Re)

    f_Im = mfem.jit.scalar(params=params)(f_exact_Im)
    g_Im = mfem.jit.scalar(params=params)(g_exact_Im)
    grad_Im = mfem.jit.vector(sdim=sdim, params=params)(grad_exact_Im)

    # ParLinearForm *b_Re(new ParLinearForm);
    b = mfem.ParComplexLinearForm(fespace, mfem.ComplexOperator.HERMITIAN)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(f_Re),
                          mfem.DomainLFIntegrator(f_Im))

    if sol >= 0:
        # if exact solution exists. Otherwise use homogeneous impedence
        # (gradp . n + i omega p = 0)
        b.AddBoundaryIntegrator(mfem.BoundaryNormalLFIntegrator(grad_Re),
                                mfem.BoundaryNormalLFIntegrator(grad_Im))
        b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(g_Re),
                                mfem.BoundaryLFIntegrator(g_Im))
    b.real().Assign(0.0)
    b.imag().Assign(0.0)

    b.Assemble()

    #  7. Set up the bilinear form (Real and Imaginary part)
    one = mfem.ConstantCoefficient(1.0)
    zero = mfem.ConstantCoefficient(0.0)
    neg_omega = mfem.ConstantCoefficient(-omega**2)

    a = mfem.ParSesquilinearForm(fespace, mfem.ComplexOperator.HERMITIAN)
    impedance = mfem.ConstantCoefficient(omega)

    bdr_attr = mfem.intArray([1]*pmesh.bdr_attributes.Max())
    imp_rest = mfem.RestrictedCoefficient(impedance, bdr_attr)

    a.AddDomainIntegrator(mfem.DiffusionIntegrator(one), None)

    # Just putting imaginary integrator with zero just to keep sparsity pattern of real
    # and imaginary part the same
    a.AddDomainIntegrator(mfem.MassIntegrator(
        neg_omega), mfem.MassIntegrator(zero))
    a.AddBoundaryIntegrator(mfem.BoundaryMassIntegrator(
        zero), mfem.BoundaryMassIntegrator(imp_rest))
    a.Assemble(0)

    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray([0] * pmesh.bdr_attributes.Max())
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    #  Solution grid function
    p_gf = mfem.ParComplexGridFunction(fespace)
    p_gf_ex = mfem.ParComplexGridFunction(fespace)

    p_Re = mfem.jit.scalar(params=params)(p_exact_Re)
    p_Im = mfem.jit.scalar(params=params)(p_exact_Im)

    p_gf.Assign(0.0)
    p_gf_ex.ProjectCoefficient(p_Re, p_Im)
    if sol >= 0:
        p_gf.ProjectBdrCoefficient(p_Re, p_Im, ess_bdr)

    Ah = mfem.OperatorHandle()
    B = mfem.Vector()
    X = mfem.Vector()

    a.FormLinearSystem(ess_tdof_list, p_gf, b, Ah, X, B)

    AZ = Ah.AsComplexHypreParMatrix()

    # ComplexMGSolver * M = new ComplexMGSolver(AZ,P,fespaces);
    # M->SetTheta(0.25);

    M = mfem.ComplexSchwarzSmoother(fespace.GetParMesh(), 0, fespace, AZ)

    gmres = mfem.GMRESSolver(MPI.COMM_WORLD)
    gmres.SetPrintLevel(1)
    gmres.SetMaxIter(2000)
    gmres.SetKDim(200)
    gmres.SetRelTol(1e-12)
    gmres.SetAbsTol(0.0)
    gmres.SetOperator(AZ)
    gmres.SetPreconditioner(M)
    gmres.Mult(B, X)

    a.RecoverFEMSolution(X, B, p_gf)

    if sol >= 0:
        h1_norm_type = 1

        L2err_Re = p_gf.real().ComputeL2Error(p_Re)
        L2err_Im = p_gf.imag().ComputeL2Error(p_Im)

        loc_H1err_Re = p_gf.real().ComputeH1Error(
            p_Re, grad_Re, one, 1.0, h1_norm_type)
        loc_H1err_Im = p_gf.imag().ComputeH1Error(
            p_Im, grad_Im, one, 1.0, h1_norm_type)

        H1err_Re = mfem.GlobalLpNorm(2.0, loc_H1err_Re, MPI.COMM_WORLD)
        H1err_Im = mfem.GlobalLpNorm(2.0, loc_H1err_Im, MPI.COMM_WORLD)

        L2error = sqrt(L2err_Re*L2err_Re + L2err_Im*L2err_Im)
        H1error = sqrt(H1err_Re*H1err_Re + H1err_Im*H1err_Im)
        if myid == 0:
            print(" || p_h - p ||_{H^1} = " + "{:g}".format(H1error))
            print(" || p_h - p ||_{L^1} = " + "{:g}".format(L2error))

    if visualization:
        keys = "keys macF\n" if dim == 3 else "keys amrRljcUUuu\n"
        if dim == 2:
            keys = "keys mrRljc\n"
        else:
            keys = "keys mc\n"

        sol_sock_re = mfem.socketstream("localhost", 19916)
        sol_sock_re.precision(8)
        sol_sock_re << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_re << "solution\n" << pmesh << p_gf.real() << keys
        sol_sock_re << "window_title 'Numerical Pressure (real part)' "
        sol_sock_re.flush()

        sol_sock_im = mfem.socketstream("localhost", 19916)
        sol_sock_im.precision(8)
        sol_sock_im << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_im << "solution\n" << pmesh << p_gf.imag() << keys
        sol_sock_re << "window_title 'Numerical Pressure (imag part)' "
        sol_sock_im.flush()


def get_helmholtz_solution_Re(x, dp, p_d2p):
    if sol == 0:   # polynomial
        if dim == 3:
            p = x[0]*(1.0 - x[0]) * x[1]*(1.0 - x[1]) * x[2]*(1.0 - x[2])
            dp[0] = (1.0 - 2.0 * x[0]) * x[1]*(1.0 - x[1]) * x[2]*(1.0 - x[2])
            dp[1] = (1.0 - 2.0 * x[1]) * x[0]*(1.0 - x[0]) * x[2]*(1.0 - x[2])
            dp[2] = (1.0 - 2.0 * x[2]) * x[0]*(1.0 - x[0]) * x[1]*(1.0 - x[1])
            d2p = (-2.0*(-1.0 + x[0]) * x[0] * (-1.0 + x[1]) * x[1]
                   - 2.0*(-1.0 + x[0]) * x[0] * (-1.0 + x[2]) * x[2]
                   - 2.0*(-1.0 + x[1]) * x[1] * (-1.0 + x[2]) * x[2])
        else:
            p = x[1] * (1.0 - x[1]) * x[0] * (1.0 - x[0])
            dp[0] = (1.0 - 2.0 * x[0]) * x[1]*(1.0 - x[1])
            dp[1] = (1.0 - 2.0 * x[1]) * x[0]*(1.0 - x[0])
            d2p = (- 2.0 * x[1] * (1.0 - x[1])
                   - 2.0 * x[0] * (1.0 - x[0]))

    elif sol == 1:  # Plane wave
        if dim == 2:
            alpha = omega/sqrt(2)
            p = cos(alpha * (x[0] + x[1]))
            dp[0] = -alpha * sin(alpha * (x[0] + x[1]))
            dp[1] = dp[0]
            d2p = -2.0 * alpha * alpha * p
        else:
            alpha = omega/sqrt(3)
            p = cos(alpha * (x[0] + x[1] + x[2]))
            dp[0] = -alpha * sin(alpha * (x[0] + x[1] + x[2]))
            dp[1] = dp[0]
            dp[2] = dp[0]
            d2p = -3.0 * alpha * alpha * p

    elif sol == 2:
        if dim == 2:
            # shift to avoid singularity
            shift = 0.1
            x0 = x[0] + shift
            x1 = x[1] + shift

            r = sqrt(x0 * x0 + x1 * x1)

            p = cos(omega * r)

            r_x = x0 / r
            r_y = x1 / r
            r_xx = (1.0 / r) * (1.0 - r_x * r_x)
            r_yy = (1.0 / r) * (1.0 - r_y * r_y)

            dp[0] = - omega * sin(omega * r) * r_x
            dp[1] = - omega * sin(omega * r) * r_y

            d2p = (-omega*omega * cos(omega * r)*r_x * r_x - omega * sin(omega*r) * r_xx
                   - omega*omega * cos(omega * r)*r_y * r_y - omega * sin(omega*r) * r_yy)
        else:
            # shift to avoid singularity
            shift = 0.1
            x0 = x[0] + shift
            x1 = x[1] + shift
            x2 = x[2] + shift

            r = sqrt(x0 * x0 + x1 * x1 + x2 * x2)

            p = cos(omega * r)

            r_x = x0 / r
            r_y = x1 / r
            r_z = x2 / r
            r_xx = (1.0 / r) * (1.0 - r_x * r_x)
            r_yy = (1.0 / r) * (1.0 - r_y * r_y)
            r_zz = (1.0 / r) * (1.0 - r_z * r_z)

            dp[0] = - omega * sin(omega * r) * r_x
            dp[1] = - omega * sin(omega * r) * r_y
            dp[2] = - omega * sin(omega * r) * r_z

            d2p = (-omega*omega * cos(omega * r)*r_x * r_x - omega * sin(omega*r) * r_xx
                   - omega*omega * cos(omega * r)*r_y *
                   r_y - omega * sin(omega*r) * r_yy
                   - omega*omega * cos(omega * r)*r_z * r_z - omega * sin(omega*r) * r_zz)

    p_d2p[0] = p
    p_d2p[1] = d2p


def get_helmholtz_solution_Im(x, dp, p_d2p):
    if sol == 0:  # polynomial
        if dim == 3:
            p = x[0]*(1.0 - x[0]) * x[1]*(1.0 - x[1]) * x[2]*(1.0 - x[2])
            dp[0] = (1.0 - 2.0 * x[0]) * x[1]*(1.0 - x[1]) * x[2]*(1.0 - x[2])
            dp[1] = (1.0 - 2.0 * x[1]) * x[0]*(1.0 - x[0]) * x[2]*(1.0 - x[2])
            dp[2] = (1.0 - 2.0 * x[2]) * x[0]*(1.0 - x[0]) * x[1]*(1.0 - x[1])
            d2p = (-2.0*(-1.0 + x[0]) * x[0] * (-1.0 + x[1]) * x[1]
                   - 2.0*(-1.0 + x[0]) * x[0] * (-1.0 + x[2]) * x[2]
                   - 2.0*(-1.0 + x[1]) * x[1] * (-1.0 + x[2]) * x[2])
        else:
            p = x[1] * (1.0 - x[1]) * x[0] * (1.0 - x[0])
            dp[0] = (1.0 - 2.0 * x[0]) * x[1]*(1.0 - x[1])
            dp[1] = (1.0 - 2.0 * x[1]) * x[0]*(1.0 - x[0])
            d2p = (- 2.0 * x[1] * (1.0 - x[1])
                   - 2.0 * x[0] * (1.0 - x[0]))
    elif sol == 1:  # plane wave
        if dim == 2:
            alpha = omega/sqrt(2)
            p = -sin(alpha * (x[0] + x[1]))
            dp[0] = -alpha * cos(alpha * (x[0] + x[1]))
            dp[1] = dp[0]
            d2p = -2.0 * alpha * alpha * p
        else:
            alpha = omega/sqrt(3)
            p = -sin(alpha * (x[0] + x[1] + x[2]))
            dp[0] = -alpha * cos(alpha * (x[0] + x[1] + x[2]))
            dp[1] = dp[0]
            dp[2] = dp[0]
            d2p = -3.0 * alpha * alpha * p
    p_d2p[0] = p
    p_d2p[1] = d2p


def p_exact_Re(x):
    p_d2p = np.empty(2)
    dp = np.empty(3)
    get_helmholtz_solution_Re(x, dp, p_d2p)
    return p_d2p[0]


def p_exact_Im(x):
    p_d2p = np.empty(2)
    dp = np.empty(3)
    get_helmholtz_solution_Im(x, dp, p_d2p)
    return p_d2p[0]


def f_exact_Re(x):
    p_d2p_re = np.empty(2)
    dp_re = np.empty(3)
    get_helmholtz_solution_Re(x, dp_re, p_d2p_re)
    p_re = p_d2p_re[0]
    d2p_re = p_d2p_re[1]
    f_re = -d2p_re - omega * omega * p_re

    return f_re


def f_exact_Im(x):
    p_d2p_im = np.empty(2)
    dp_im = np.empty(3)

    get_helmholtz_solution_Im(x, dp_im, p_d2p_im)
    p_im = p_d2p_im[0]
    d2p_im = p_d2p_im[1]

    f_im = -d2p_im - omega * omega * p_im
    return f_im


def grad_exact_Re(x, ret):
    p_d2p = np.empty(2)
    dp = np.empty(3)
    get_helmholtz_solution_Re(x, dp, p_d2p)
    for i in range(dim):
        ret[i] = dp[i]


def grad_exact_Im(x, ret):
    p_d2p = np.empty(2)
    dp = np.empty(3)
    get_helmholtz_solution_Im(x, dp, p_d2p)
    for i in range(dim):
        ret[i] = dp[i]

# define impedence coefficient: i omega p


def g_exact_Re(x):
    p_d2p = np.empty(2)
    dp = np.empty(3)
    get_helmholtz_solution_Im(x, dp, p_d2p)
    return -omega * p_d2p[0]


def g_exact_Im(x):
    p_d2p = np.empty(2)
    dp = np.empty(3)
    get_helmholtz_solution_Re(x, dp, p_d2p)
    return omega * p_d2p[0]


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='complex MG')

    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Elemet order")
    parser.add_argument('-nd', '--ndim',
                        action='store', default=2, type=int,
                        help='Problem space dimension')
    parser.add_argument("-sol", "--exact_sol",
                        action='store', default=1, type=int,
                        help="Exact solution flag - 0:polynomial, " +
                        "1: plane wave, -1: unknown exact")
    parser.add_argument("-k", "--wavelengths",
                        action='store', default=0.5, type=float,
                        help="Number of wavelengths.")
    parser.add_argument('-length', '--length',
                        action='store', default=1, type=float,
                        help='length of domain')
    parser.add_argument('-href', '--href',
                        action='store', default=1, type=int,
                        help="Number of geometric refinements done " +
                        "prior to order refinements.")
    parser.add_argument('-initref', '--initref',
                        action='store', default=1, type=int,
                        help='initial refinement')
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    mesh_file = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    run(order=args.order,
        href=args.href,
        mesh_file=mesh_file,
        nd=args.ndim,
        initref=args.initref,
        length=args.length,
        k=args.wavelengths,
        sol=args.exact_sol,
        visualization=args.visualization)
