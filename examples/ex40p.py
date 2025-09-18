'''
   MFEM example 40p (converted from ex40p.cpp)

   See c++ version in the MFEM library for more detail

   Sample runs: mpirun -np  4 python ex40.py -step 10.0 -gr 2.0
                mpirun -np  4 python ex40.py -step 10.0 -gr 2.0 -o 3 -r 1
                mpirun -np  4 python ex40.py -step 10.0 -gr 2.0 -r 4 -m ../data/l-shape.mesh
                mpirun -np  4 python ex40.py -step 10.0 -gr 2.0 -r 2 -m ../data/fichera.mesh

   Description: This example code demonstrates how to use MFEM to solve the
                eikonal equation,

                        |∇𝑢| = 1 in Ω,  𝑢 = 0 on ∂Ω.

                The viscosity solution of this problem coincides with the unique optimum
                of the nonlinear program

                     maximize ∫_Ω 𝑢 d𝑥 subject to |∇𝑢| ≤ 1 in Ω, 𝑢 = 0 on ∂Ω,    (⋆)

                which is the foundation for method implemented below.

                Following the proximal Galerkin methodology [1,2] (see also Example
                36), we construct a Legendre function for the closed unit ball
                𝐵₁ := {𝑥 ∈ Rⁿ | |𝑥| ≤ 1}. Our choice is the Hellinger entropy,

                      R(𝑥) = −( 1 − |𝑥|² )^{1/2},

                although other choices are possible, each leading to a slightly
                different algorithm. We then adaptively regularize the optimization
                problem (⋆) with the Bregman divergence of the Hellinger entropy,

                   maximize  ∫_Ω 𝑢 d𝑥 - αₖ⁻¹ D(∇𝑢,∇𝑢ₖ₋₁)  subject to  𝑢 = 0 on Ω.

                This results in a sequence of functions ( 𝜓ₖ , 𝑢ₖ ),

                        𝑢ₖ → 𝑢,    𝜓ₖ/|𝜓ₖ| → ∇𝑢    as k → ∞,

                defined by the nonlinear saddle-point problems

                 Find 𝜓ₖ ∈ H(div,Ω) and 𝑢ₖ ∈ L²(Ω) such that
                 ( (∇R)⁻¹(𝜓ₖ) , τ ) + ( 𝑢ₖ , ∇⋅τ ) = 0                     ∀ τ ∈ H(div,Ω)
                 ( ∇⋅𝜓ₖ , v )                     = ( ∇⋅𝜓ₖ₋₁ - αₖ , v )    ∀ v ∈ L²(Ω)

                where (∇R)⁻¹(𝜓) = 𝜓 / ( 1 + |𝜓|² )^{1/2} and αₖ = α₀rᵏ, where r ≥ 1
                is a prescribed growth rate. (r = 1 is the most stable.) The
                saddle-point problems are solved using a damped quasi-Newton method
                with a tunable regularization parameter 0 ≤ ϵ << 1.

                [1] Keith, B. and Surowiec, T. (2024) Proximal Galerkin: A structure-
                    preserving finite element method for pointwise bound constraints.
                    Foundations of Computational Mathematics, 1–97.
                [2] Dokken, J., Farrell, P., Keith, B., Papadopoulos, I., and
                    Surowiec, T. (2025) The latent variable proximal point algorithm
                    for variational problems with inequality constraints. (To appear.)

'''
import os
from os.path import expanduser, join
import numpy as np

import mfem.par as mfem

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(meshfile="",
        order=1,
        max_it=5,
        ref_levels=3,
        alpha=1.0,
        growth_rate=1.0,
        newton_scaling=0.8,
        eps=1e-6,
        tol=1e-4,
        max_alpha=1e2,
        max_psi=1e2,
        eps2=1e-1,
        visualization=True):

    # 2. Read the mesh from the mesh file.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    # 3. Postprocess the mesh.
    # 3A. Refine the mesh to increase the resolution.
    for i in range(ref_levels):
        mesh.UniformRefinement()

    # 3B. Interpolate the geometry after refinement to control geometry error.
    # NOTE: Minimum second-order interpolation is used to improve the accuracy.
    curvature_order = max(order, 2)
    mesh.SetCurvature(curvature_order)

    # 3C. Compute the maximum mesh size.
    hmax = 0.0
    for i in range(mesh.GetNE()):
        hmax = max(mesh.GetElementSize(i, 1), hmax)

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh

    # 4. Define the necessary finite element spaces on the mesh.
    RTfec = mfem.RT_FECollection(order, dim)
    RTfes = mfem.ParFiniteElementSpace(pmesh, RTfec)

    L2fec = mfem.L2_FECollection(order, dim)
    L2fes = mfem.ParFiniteElementSpace(pmesh, L2fec)

    dofs_rt = RTfes.GlobalTrueVSize()
    dofs_l2 = L2fes.GlobalTrueVSize()
    if myid == 0:
        print("Number of H(div) dofs: " + str(dofs_rt))
        print("Number of L² dofs: " + str(dofs_l2))

    # 5. Define the offsets for the block matrices
    offsets = mfem.intArray([0, RTfes.GetVSize(), L2fes.GetVSize()])
    offsets.PartialSum()
    toffsets = mfem.intArray([0, RTfes.GetTrueVSize(), L2fes.GetTrueVSize()])
    toffsets.PartialSum()

    x = mfem.BlockVector(offsets)
    rhs = mfem.BlockVector(offsets)
    x.Assign(0.0)
    rhs.Assign(0.0)
    tx = mfem.BlockVector(toffsets)
    trhs = mfem.BlockVector(toffsets)
    tx.Assign(0.0)
    trhs.Assign(0.0)

    # 6. Define the solution vectors as a finite element grid functions
    #    corresponding to the fespaces.

    u_gf = mfem.ParGridFunction()
    delta_psi_gf = mfem.ParGridFunction()
    delta_psi_gf.MakeRef(RTfes, x, offsets[0])
    u_gf.MakeRef(L2fes, x, offsets[1])

    psi_old_gf = mfem.ParGridFunction(RTfes)
    psi_gf = mfem.ParGridFunction(RTfes)
    u_old_gf = mfem.ParGridFunction(L2fes)

    # 7. Define initial guesses for the solution variables.
    delta_psi_gf.Assign(0.0)
    psi_gf.Assign(0.0)
    u_gf.Assign(0.0)
    psi_old_gf.Assign(0.0)
    u_old_gf.Assign(0.0)

    # 8. Prepare for glvis output.
    # 14. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)

    # 9. Coefficients to be used later.
    zero_cf = mfem.ConstantCoefficient(0.0)
    one_cf = mfem.ConstantCoefficient(1.0)
    zero_vec = mfem.Vector([0.0]*sdim)
    zero_vec_cf = mfem.VectorConstantCoefficient(zero_vec)
    neg_alpha_cf = mfem.ConstantCoefficient(-1.0*alpha)

    psigf_cf = mfem.VectorGridFunctionCoefficient(psi_gf)

    @mfem.jit.vector(shape=(sdim,), dependency=(psigf_cf,))
    def Z(_ptx, psi_vals):
        norm = np.linalg.norm(psi_vals)
        phi = 1.0 / np.sqrt(1.0 + norm*norm)
        vvec = psi_vals*phi
        return vvec

    @mfem.jit.matrix(shape=(sdim, sdim), dependency=(psigf_cf,))
    def DZ(_ptx, psi_vals):
        norm = np.linalg.norm(psi_vals)
        phi = 1.0 / np.sqrt(1.0 + norm*norm)

        kmat = np.zeros((sdim, sdim))
        for i in range(sdim):
            kmat[i, i] = phi + eps
            for j in range(sdim):
                kmat[i, j] -= psi_vals[i]*psi_vals[j] * phi**3
        return kmat

    neg_Z = mfem.ScalarVectorProductCoefficient(-1.0, Z)
    div_psi_cf = mfem.DivergenceGridFunctionCoefficient(psi_gf)
    div_psi_old_cf = mfem.DivergenceGridFunctionCoefficient(psi_old_gf)
    psi_old_minus_psi = mfem.SumCoefficient(
        div_psi_old_cf, div_psi_cf, 1.0, -1.0)

    # 10. Assemble constant matrices/vectors to avoid reassembly in the loop.
    b0 = mfem.ParLinearForm()
    b1 = mfem.ParLinearForm()
    b0.MakeRef(RTfes, rhs.GetBlock(0), 0)
    b1.MakeRef(L2fes, rhs.GetBlock(1), 0)

    b0.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(neg_Z))
    b1.AddDomainIntegrator(mfem.DomainLFIntegrator(neg_alpha_cf))
    b1.AddDomainIntegrator(mfem.DomainLFIntegrator(psi_old_minus_psi))

    a00 = mfem.ParBilinearForm(RTfes)
    a00.AddDomainIntegrator(mfem.VectorFEMassIntegrator(DZ))

    a10 = mfem.ParMixedBilinearForm(RTfes, L2fes)
    a10.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator())
    a10.Assemble()
    a10.Finalize()
    A10 = a10.ParallelAssemble()
    A01 = A10.Transpose()

    vol_form = mfem.ParLinearForm(L2fes)
    one_gf = mfem.ParGridFunction(L2fes)
    one_gf.Assign(1.0)
    vol_form.AddDomainIntegrator(mfem.DomainLFIntegrator(one_cf))
    vol_form.Assemble()
    domain_volume = vol_form(one_gf)

    # 11. Iterate.
    total_iterations = 0
    increment_u = 0.1
    u_tmp = mfem.ParGridFunction(L2fes)

    for k in range(max_it):
        u_tmp.Assign(u_old_gf)

        if myid == 0:
            print("\nOUTER ITERATION " + str(k+1))

        for j in range(5):
            total_iterations += 1

            b0.Assemble()
            b0.ParallelAssemble(trhs.GetBlock(0))

            b1.Assemble()
            b1.ParallelAssemble(trhs.GetBlock(1))

            a00.Assemble(0)
            a00.Finalize(0)
            A00 = a00.ParallelAssemble()

            # Construct Schur-complement preconditioner
            A00_diag = mfem.HypreParVector(MPI.COMM_WORLD,
                                           A00.GetGlobalNumRows(),
                                           A00.GetRowStarts())
            A00.GetDiag(A00_diag)
            S_tmp = mfem.HypreParMatrix(A01)
            S_tmp.InvScaleRows(A00_diag)
            A00_diag.Reciprocal()
            S = mfem.ParMult(A10, S_tmp, True)

            # Python note:
            #    owns_blocks should be 0 in Python
            #    because wrapper class will delete blocks
            prec = mfem.BlockDiagonalPreconditioner(toffsets)
            prec.owns_blocks = 0

            P00 = mfem.HypreBoomerAMG(A00)
            P00.SetPrintLevel(0)
            P11 = mfem.HypreBoomerAMG(S)
            P11.SetPrintLevel(0)
            prec.SetDiagonalBlock(0, P00)
            prec.SetDiagonalBlock(1, P11)

            A = mfem.BlockOperator(toffsets)
            A.SetBlock(0, 0, A00)
            A.SetBlock(1, 0, A10)
            A.SetBlock(0, 1, A01)

            minres = mfem.MINRESSolver(MPI.COMM_WORLD)
            minres.SetPrintLevel(-1)
            minres.SetRelTol(1e-12)
            minres.SetMaxIter(10000)
            minres.SetOperator(A)
            minres.SetPreconditioner(prec)
            minres.Mult(trhs, tx)

            del S
            del A00

            delta_psi_gf.SetFromTrueDofs(tx.GetBlock(0))
            u_gf.SetFromTrueDofs(tx.GetBlock(1))

            u_tmp -= u_gf  # u_tmp = u_tmp - u_gf

            Newton_update_size = u_tmp.ComputeL2Error(zero_cf)
            u_tmp.Assign(u_gf)

            # Damped Newton update
            psi_gf.Add(newton_scaling, delta_psi_gf)
            a00.Update()

            if visualization:
                sol_sock << "parallel " << num_procs << " " << myid << "\n"
                sol_sock << "solution\n" << pmesh << u_gf
                sol_sock << "window_title 'Discrete solution'"
                sol_sock.flush()

            if myid == 0:
                print("Newton_update_size = " +
                      "{:g}".format(Newton_update_size))

            if newton_scaling*Newton_update_size < increment_u:
                break

        u_tmp.Assign(u_gf)
        u_tmp -= u_old_gf
        increment_u = u_tmp.ComputeL2Error(zero_cf)

        if myid == 0:
            print("Number of Newton iterations = " + str(j+1))
            print("Increment (|| uₕ - uₕ_prvs||) = " +
                  "{:g}".format(increment_u))

        u_old_gf.Assign(u_gf)
        psi_old_gf.Assign(psi_gf)
        alpha *= max(growth_rate, 1.0)

        #  Safeguard 1: Stop alpha from growing too large
        alpha = min(alpha, max_alpha)

        # Safeguard 2: Stop |ψ| from growing too large
        norm_psi = psi_old_gf.ComputeL1Error(zero_vec_cf)/domain_volume
        if norm_psi > max_psi:
            # Additional entropy regularization
            neg_alpha_cf.constant = -alpha/(1.0 + eps2 * alpha * hmax)
            psi_old_minus_psi.SetAlpha(1.0/(1.0 + eps2 * alpha * hmax))
        else:
            neg_alpha_cf.constant = -alpha

        if increment_u < tol or k == max_it-1:
            break

        alpha *= max(growth_rate, 1.0)
        neg_alpha_cf.constant = -alpha

    if myid == 0:
        print("\n Outer iterations: " + str(k+1) +
              "\n Total iterations: " + str(total_iterations) +
              "\n Total dofs:       " + str(RTfes.GetTrueVSize() + L2fes.GetTrueVSize()))


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex40 (Eikonal queation)')
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).")
    parser.add_argument('-r', '--refs',
                        action='store', default=3, type=int,
                        help="Number of h-refinements.")
    parser.add_argument('-mi', '--max-it',
                        action='store', default=5, type=int,
                        help="Maximum number of iterations")
    parser.add_argument("-tol", "--tol",
                        action='store', default=1e-4, type=float,
                        help="Stopping criteria based on the difference between.")
    parser.add_argument('-step', '--step',
                        action='store', default=1.0, type=float,
                        help="Initial size alpha")
    parser.add_argument("-gr", "--growth-rate",
                        action='store', default=1.0, type=float,
                        help="Growth rate of the step size alpha")
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        default=False,
                        help='Disable or disable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = not args.no_visualization
    in_alpha = args.step

    run(meshfile=meshfile,
        order=args.order,
        max_it=args.max_it,
        ref_levels=args.refs,
        alpha=in_alpha,
        growth_rate=args.growth_rate,
        newton_scaling=0.8,
        eps=1e-6,
        tol=args.tol,
        visualization=visualization)
