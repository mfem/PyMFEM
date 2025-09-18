'''
   MFEM example 40 (converted from ex40.cpp)

   See c++ version in the MFEM library for more detail

   Sample runs: python ex40.py -step 10.0 -gr 2.0
                python ex40.py -step 10.0 -gr 2.0 -o 3 -r 1
                python ex40.py -step 10.0 -gr 2.0 -r 4 -m ../data/l-shape.mesh
                python ex40.py -step 10.0 -gr 2.0 -r 2 -m ../data/fichera.mesh

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

import mfem.ser as mfem

if hasattr(mfem, "UMFPackSolver"):
    use_umfpack = True if not args.no_use_umfpack else False
else:
    use_umfpack = False


def run(meshfile="",
        order=1,
        max_it=5,
        ref_levels=3,
        alpha=1.0,
        growth_rate=1.0,
        newton_scaling=0.8,
        eps=1e-6,
        tol=1e-4,
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

    # 4. Define the necessary finite element spaces on the mesh.
    RTfec = mfem.RT_FECollection(order, dim)
    RTfes = mfem.FiniteElementSpace(mesh, RTfec)

    L2fec = mfem.L2_FECollection(order, dim)
    L2fes = mfem.FiniteElementSpace(mesh, L2fec)

    print("Number of H(div) dofs: " + str(RTfes.GetTrueVSize()))
    print("Number of L² dofs: " + str(L2fes.GetTrueVSize()))

    # 5. Define the offsets for the block matrices
    offsets = mfem.intArray([0, RTfes.GetVSize(), L2fes.GetVSize()])
    offsets.PartialSum()

    x = mfem.BlockVector(offsets)
    rhs = mfem.BlockVector(offsets)
    x.Assign(0.0)
    rhs.Assign(0.0)

    # 6. Define the solution vectors as a finite element grid functions
    #    corresponding to the fespaces.

    u_gf = mfem.GridFunction()
    delta_psi_gf = mfem.GridFunction()
    delta_psi_gf.MakeRef(RTfes, x, offsets[0])
    u_gf.MakeRef(L2fes, x, offsets[1])

    psi_old_gf = mfem.GridFunction(RTfes)
    psi_gf = mfem.GridFunction(RTfes)
    u_old_gf = mfem.GridFunction(L2fes)

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
    neg_alpha_cf = mfem.ConstantCoefficient(-1.0*alpha)
    zero_cf = mfem.ConstantCoefficient(0.0)

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
    b0 = mfem.LinearForm()
    b1 = mfem.LinearForm()
    b0.MakeRef(RTfes, rhs.GetBlock(0), 0)
    b1.MakeRef(L2fes, rhs.GetBlock(1), 0)

    b0.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(neg_Z))
    b1.AddDomainIntegrator(mfem.DomainLFIntegrator(neg_alpha_cf))
    b1.AddDomainIntegrator(mfem.DomainLFIntegrator(psi_old_minus_psi))

    a00 = mfem.BilinearForm(RTfes)
    a00.AddDomainIntegrator(mfem.VectorFEMassIntegrator(DZ))

    a10 = mfem.MixedBilinearForm(RTfes, L2fes)
    a10.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator())
    a10.Assemble()
    a10.Finalize()
    A10 = a10.SpMat()
    A01 = mfem.Transpose(A10)

    # 11. Iterate.
    total_iterations = 0
    increment_u = 0.1
    u_tmp = mfem.GridFunction(L2fes)

    for k in range(max_it):
        u_tmp.Assign(u_old_gf)

        print("\nOUTER ITERATION " + str(k+1))

        for j in range(5):
            total_iterations += 1

            b0.Assemble()
            b1.Assemble()

            a00.Assemble(0)
            a00.Finalize(0)
            A00 = a00.SpMat()

            # Construct Schur-complement preconditioner
            A00_diag = mfem.Vector(a00.Height())
            A00.GetDiag(A00_diag)
            A00_diag.Reciprocal()
            S = mfem.Mult_AtDA(A01, A00_diag)

            # Python note:
            #    owns_blocks should be 0 in Python
            #    because wrapper class will delete blocks
            prec = mfem.BlockDiagonalPreconditioner(offsets)
            prec.owns_blocks = 0

            prec.SetDiagonalBlock(0, mfem.DSmoother(A00))

            if not use_umfpack:
                prec.SetDiagonalBlock(1, mfem.GSSmoother(S))
            else:
                prec.SetDiagonalBlock(1, mfem.UMFPackSolver(S))

            A = mfem.BlockOperator(offsets)
            A.SetBlock(0, 0, A00)
            A.SetBlock(1, 0, A10)
            A.SetBlock(0, 1, A01)

            mfem.MINRES(A, prec, rhs, x, 0, 2000, 1e-12)

            del S

            u_tmp -= u_gf  # u_tmp = u_tmp - u_gf

            Newton_update_size = u_tmp.ComputeL2Error(zero_cf)
            u_tmp.Assign(u_gf)

            # Damped Newton update
            psi_gf.Add(newton_scaling, delta_psi_gf)
            a00.Update()

            if visualization:
                sol_sock << "solution\n" << mesh << u_gf
                sol_sock << "window_title 'Discrete solution'"
                sol_sock.flush()

            print("Newton_update_size = " + "{:g}".format(Newton_update_size))

            if Newton_update_size < increment_u:
                break

        u_tmp.Assign(u_gf)
        u_tmp -= u_old_gf
        increment_u = u_tmp.ComputeL2Error(zero_cf)

        print("Number of Newton iterations = " + str(j+1))
        print("Increment (|| uₕ - uₕ_prvs||) = " + "{:g}".format(increment_u))

        u_old_gf.Assign(u_gf)
        psi_old_gf.Assign(psi_gf)

        if increment_u < tol or k == max_it-1:
            break

        alpha *= max(growth_rate, 1.0)
        neg_alpha_cf.constant = -alpha

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
