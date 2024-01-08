'''
   PyMFEM example 36p

   See c++ version in the MFEM library for more detail

   Sample runs:  python ex36.py -o 2
                 python ex36.py -o 2 -r 4


'''
import os
from os.path import expanduser, join

from numpy import sqrt, log, exp
import numpy as np
from numba import njit
from numba.types import float64


import mfem.par as mfem
from mfem.par import intArray, doubleArray
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)

visualization = True


def run(refs=3,
        order=1,
        max_it=10,
        tol=1e-5,
        alpha=1.0):

    # 2. Read the mesh from the mesh file.
    mesh_file = "../data/disc-nurbs.mesh"
    mesh_file = expanduser(
        join(os.path.dirname(__file__), *mesh_file.split("/")))
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()

    # 3. Postprocess the mesh.
    # 3A. Refine the mesh to increase the resolution.
    for l in range(refs):
        mesh.UniformRefinement()

    # 3B. Interpolate the geometry after refinement to control geometry error.
    # NOTE: Minimum second-order interpolation is used to improve the accuracy.
    curvature_order = max((order, 2))
    mesh.SetCurvature(curvature_order)

    # 3C. Rescale the domain to a unit circle (radius = 1).
    nodes = mesh.GetNodes()
    scale = 2*sqrt(2)
    nodes /= scale

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    mesh.Clear()

    # 4. Define the necessary finite element spaces on the mesh.
    H1fec = mfem.H1_FECollection(order+1, dim)
    H1fes = mfem.ParFiniteElementSpace(pmesh, H1fec)

    L2fec = mfem.L2_FECollection(order-1, dim)
    L2fes = mfem.ParFiniteElementSpace(pmesh, L2fec)

    num_dofs_H1 = H1fes.GlobalTrueVSize()
    num_dofs_L2 = L2fes.GlobalTrueVSize()
    if myid == 0:
        print("Number of H1 finite element unknowns: " + str(num_dofs_H1))
        print("Number of L2 finite element unknowns: " + str(num_dofs_L2))

    offsets = mfem.intArray((0, H1fes.GetVSize(), L2fes.GetVSize()))
    offsets.PartialSum()
    toffsets = mfem.intArray((0, H1fes.GetTrueVSize(), L2fes.GetTrueVSize()))
    toffsets.PartialSum()

    x = mfem.BlockVector(offsets)
    rhs = mfem.BlockVector(offsets)
    x.Assign(0.0)
    rhs.Assign(0.0)
    tx = mfem.BlockVector(toffsets)
    trhs = mfem.BlockVector(toffsets)
    tx.Assign(0.0)
    trhs.Assign(0.0)

    # 5. Determine the list of true (i.e. conforming) essential boundary dofs.
    empty = mfem.intArray()
    ess_tdof_list = mfem.intArray()
    if pmesh.bdr_attributes.Size() > 0:
        ess_bdr = mfem.intArray([1]*pmesh.bdr_attributes.Max())
        H1fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 6. Define an initial guess for the solution.
    @njit(float64(float64[:]))
    def IC_func(x):
        r0 = 1.0
        rr = 0.0
        for i in range(len(x)):
            rr += x[i]*x[i]
        return r0*r0 - rr

    one = mfem.ConstantCoefficient(1.0)
    zero = mfem.ConstantCoefficient(0.0)

    # 7. Define the solution vectors as a finite element grid functions
    #    corresponding to the fespaces.
    u_gf = mfem.ParGridFunction()
    delta_psi_gf = mfem.ParGridFunction()
    u_gf.MakeRef(H1fes, x, offsets[0])
    delta_psi_gf.MakeRef(L2fes, x, offsets[1])
    delta_psi_gf.Assign(0.0)

    u_old_gf = mfem.ParGridFunction(H1fes)
    psi_old_gf = mfem.ParGridFunction(L2fes)
    psi_gf = mfem.ParGridFunction(L2fes)
    u_old_gf.Assign(0.0)
    psi_old_gf.Assign(0.0)

    # 8. Define the function coefficients for the solution and use them to
    #    initialize the initial guess
    @mfem.jit.scalar
    def exact_coef(ptx):
        return exact_solution_obstacle(ptx[0], ptx[1])

    @mfem.jit.vector(shape=(dim, ))
    def exact_grad_coef(ptx):
        return exact_solution_gradient_obstacle(ptx[0], ptx[1])

    @mfem.jit.scalar
    def IC_coef(ptx):
        return IC_func(ptx)
    f = mfem.ConstantCoefficient(0.0)

    @mfem.jit.scalar
    def obstacle(ptx):
        return spherical_obstacle(ptx[0], ptx[1])

    u_gf.ProjectCoefficient(IC_coef)
    u_old_gf.GetDataArray()[:] = u_gf.GetDataArray()

    # 9. Initialize the slack variable ψₕ = ln(uₕ)
    c_u_gf = mfem.GridFunctionCoefficient(u_gf)

    @mfem.jit.scalar(dependency=(c_u_gf, obstacle))
    def ln_u(ptx, u, obstacle):
        min_val = -36.
        val = u - obstacle
        return max((min_val, log(val)))

    psi_gf.ProjectCoefficient(ln_u)
    psi_old_gf.GetDataArray()[:] = psi_gf.GetDataArray()

    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)

    # Note: we define exp_psi outside the iteration loop to avoid
    #       numba-jit during the loop.
    #       GridFunction referred from c_psi_gf is updated during
    #       the iteration.
    c_psi_gf = mfem.GridFunctionCoefficient()

    @mfem.jit.scalar(dependency=(c_psi_gf, ))
    def exp_psi(ptx, val):
        min_val = 0
        max_val = 1e36
        return min(max_val, max(min_val, exp(val)))

    # 10. Iterate
    total_iterations = 0
    increment_u = 0.1
    for k in range(max_it):
        u_tmp = mfem.ParGridFunction(H1fes)
        u_tmp.GetDataArray()[:] = u_old_gf.GetDataArray()

        if myid == 0:
            print("\nOUTER ITERATION " + str(k+1))

        for j in range(10):
            total_iterations += 1

            alpha_cf = mfem.ConstantCoefficient(alpha)

            b0 = mfem.ParLinearForm()
            b1 = mfem.ParLinearForm()
            b0.Update(H1fes, rhs.GetBlock(0), 0)
            b1.Update(L2fes, rhs.GetBlock(1), 0)

            c_psi_gf.SetGridFunction(psi_gf)
            neg_exp_psi = mfem.ProductCoefficient(-1.0, exp_psi)
            grad_u_old = mfem.GradientGridFunctionCoefficient(u_old_gf)
            alpha_f = mfem.ProductCoefficient(alpha, f)
            psi_cf = mfem.GridFunctionCoefficient(psi_gf)
            psi_old_cf = mfem.GridFunctionCoefficient(psi_old_gf)
            psi_old_minus_psi = mfem.SumCoefficient(
                psi_old_cf, psi_cf, 1.0, -1.0)

            b0.AddDomainIntegrator(mfem.DomainLFIntegrator(alpha_f))
            b0.AddDomainIntegrator(mfem.DomainLFIntegrator(psi_old_minus_psi))
            b0.Assemble()

            b1.AddDomainIntegrator(mfem.DomainLFIntegrator(exp_psi))
            b1.AddDomainIntegrator(mfem.DomainLFIntegrator(obstacle))
            b1.Assemble()

            a00 = mfem.ParBilinearForm(H1fes)
            a00.SetDiagonalPolicy(mfem.Operator.DIAG_ONE)
            a00.AddDomainIntegrator(mfem.DiffusionIntegrator(alpha_cf))
            a00.Assemble()
            A00 = mfem.HypreParMatrix()
            a00.FormLinearSystem(ess_tdof_list, x.GetBlock(0), rhs.GetBlock(0),
                                 A00, tx.GetBlock(0), trhs.GetBlock(0))

            a10 = mfem.ParMixedBilinearForm(H1fes, L2fes)
            a10.AddDomainIntegrator(mfem.MixedScalarMassIntegrator())
            a10.Assemble()
            A10 = mfem.HypreParMatrix()
            a10.FormRectangularLinearSystem(ess_tdof_list, empty, x.GetBlock(0),
                                            rhs.GetBlock(1),
                                            A10, tx.GetBlock(0), trhs.GetBlock(1))

            A01 = A10.Transpose()

            a11 = mfem.ParBilinearForm(L2fes)
            a11.AddDomainIntegrator(mfem.MassIntegrator(neg_exp_psi))
            # NOTE: Shift the spectrum of the Hessian matrix for additional
            #       stability (Quasi-Newton).
            eps_cf = mfem.ConstantCoefficient(-1e-6)
            if (order == 1):
                # NOTE: ∇ₕuₕ = 0 for constant functions.
                #       Therefore, we use the mass matrix to shift the spectrum
                a11.AddDomainIntegrator(mfem.MassIntegrator(eps_cf))
            else:
                a11.AddDomainIntegrator(mfem.DiffusionIntegrator(eps_cf))

            a11.Assemble()
            A11 = mfem.HypreParMatrix()
            a11.FormSystemMatrix(empty, A11)

            A = mfem.BlockOperator(toffsets)
            A.SetBlock(0, 0, A00)
            A.SetBlock(1, 0, A10)
            A.SetBlock(0, 1, A01)
            A.SetBlock(1, 1, A11)

            prec = mfem.BlockDiagonalPreconditioner(toffsets)
            P00 = mfem.HypreBoomerAMG(A00)
            P11 = mfem.HypreBoomerAMG(A11)
            P00.SetPrintLevel(0)
            P11.SetPrintLevel(0)
            prec.SetDiagonalBlock(0, P00)
            prec.SetDiagonalBlock(1, P11)
            prec.owns_blocks = 0

            gmres = mfem.GMRESSolver(MPI.COMM_WORLD)
            gmres.SetPrintLevel(-1)
            gmres.SetRelTol(1e-8)
            gmres.SetMaxIter(20000)
            gmres.SetKDim(500)
            gmres.SetOperator(A)
            gmres.SetPreconditioner(prec)
            gmres.Mult(trhs, tx)

            u_gf.SetFromTrueDofs(tx.GetBlock(0))
            delta_psi_gf.SetFromTrueDofs(tx.GetBlock(1))

            u_tmp -= u_gf
            Newton_update_size = u_tmp.ComputeL2Error(zero)
            u_tmp.Assign(u_gf)
            # Note: to copy GridFunction data, use either Assign or
            #       numpy representation.
            #    u_tmp.Assign(u_gf)
            #    u_tmp.GetDataArray()[:] = u_gf.GetDataArray()

            gamma = 1.0
            delta_psi_gf *= gamma
            psi_gf += delta_psi_gf

            if visualization:
                sol_sock << "parallel " << num_procs << " " << myid << "\n"
                sol_sock << "solution\n" << pmesh << u_gf
                sol_sock << "window_title 'Discrete solution'"
                sol_sock.flush()

            if myid == 0:
                print("Newton_update_size = {:g}".format(Newton_update_size))

            del A01

            if Newton_update_size < increment_u:
                break

        u_tmp.Assign(u_gf)
        u_tmp -= u_old_gf
        increment_u = u_tmp.ComputeL2Error(zero)

        if myid == 0:
            print("Number of Newton iterations = " + str(j+1))
            print("Increment (|| uₕ - uₕ_prvs||) = " +
                  "{:g}".format(increment_u))

        u_old_gf.Assign(u_gf)
        psi_old_gf.Assign(psi_gf)

        if increment_u < tol or k == max_it-1:
            break

        H1_error = u_gf.ComputeH1Error(exact_coef, exact_grad_coef)
        if myid == 0:
            print("H1-error  (|| u - uₕᵏ||)       = " +
                  "{:g}".format(H1_error))

    if myid == 0:
        print("Outer iterations: " + str(k+1))
        print(" Total iterations: " + str(total_iterations))
        print(" Total dofs:       " +
              str(H1fes.GlobalTrueVSize() + L2fes.GlobalTrueVSize()))

    # 11. Exact solution.
    if visualization:
        err_sock = mfem.socketstream("localhost", 19916)
        err_sock.precision(8)

        error_gf = mfem.ParGridFunction(H1fes)
        error_gf.ProjectCoefficient(exact_coef)
        error_gf -= u_gf

        err_sock << "parallel " << num_procs << " " << myid << "\n"
        err_sock << "solution\n" << pmesh << error_gf << "window_title 'Error'"
        err_sock.flush()

    L2_error = u_gf.ComputeL2Error(exact_coef)
    H1_error = u_gf.ComputeH1Error(exact_coef, exact_grad_coef)

    c_psi_gf = mfem.GridFunctionCoefficient(psi_gf)

    @mfem.jit.scalar(dependency=(c_psi_gf, obstacle))
    def u_alt_cf(ptx, val, obstacle):
        min_val = 0
        max_val = 1e36
        return min(max_val, max(min_val, exp(val)+obstacle))

    u_alt_gf = mfem.ParGridFunction(L2fes)
    u_alt_gf.ProjectCoefficient(u_alt_cf)
    L2_error_alt = u_alt_gf.ComputeL2Error(exact_coef)

    if myid == 0:
        print("\n Final L2-error (|| u - uₕ||)          = " +
              "{:g}".format(L2_error))
        print(" Final H1-error (|| u - uₕ||)          = " +
              "{:g}".format(H1_error))
        print(" Final L2-error (|| u - ϕ - exp(ψₕ)||) = " +
              "{:g}".format(L2_error_alt))


@njit(float64(float64, float64))
def spherical_obstacle(x, y):
    r = sqrt(x*x + y*y)
    r0 = 0.5
    beta = 0.9

    b = r0*beta
    tmp = sqrt(r0*r0 - b*b)
    B = tmp + b*b/tmp
    C = -b/tmp

    if r > b:
        return B + r * C
    else:
        return sqrt(r0*r0 - r*r)


@njit(float64(float64, float64))
def exact_solution_obstacle(x, y):
    r = sqrt(x*x + y*y)
    r0 = 0.5
    a = 0.348982574111686
    A = -0.340129705945858

    if r > a:
        return A * log(r)
    else:
        return sqrt(r0*r0-r*r)


@njit(float64[:](float64, float64))
def exact_solution_gradient_obstacle(x, y):
    r = sqrt(x*x + y*y)
    r0 = 0.5
    a = 0.348982574111686
    A = -0.340129705945858

    grad = np.zeros(2, dtype=np.float64)
    if r > a:
        grad[0] = A * x / (r*r)
        grad[1] = A * y / (r*r)
    else:
        grad[0] = - x / sqrt(r0*r0 - r*r)
        grad[1] = - y / sqrt(r0*r0 - r*r)
    return grad


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex36p (Obstacle problem)')

    parser.add_argument("-o", "--order",
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).")
    parser.add_argument('-r', '--refs',
                        action='store', default=3, type=int,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument('-mi', '--max-it',
                        action='store', default=10, type=int,
                        help="Maximum number of iterations")
    help = "\n".join(("Stopping criteria based on the difference between",
                      "successive solution updates"))
    parser.add_argument('-tol', '--tol',
                        action='store', default=1e-5, type=float,
                        help=help)
    parser.add_argument('-step', '--step',
                        action='store', default=1.0, type=float,
                        help="Step size alpha")
    parser.add_argument("-no-vis", "--no-visualization",
                        action='store_true', default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    globals()["visualization"] = not args.no_visualization

    run(refs=args.refs,
        order=args.order,
        max_it=args.max_it,
        tol=args.tol,
        alpha=args.step)
