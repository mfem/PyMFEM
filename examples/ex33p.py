'''
   MFEM example 33

   See c++ version in the MFEM library for more detail

   Sample runs:  mpirun -np 2 ex33p.py -m ../data/square-disc.mesh -alpha 0.33 -o 2
                 mpirun -np 2 ex33p.py -m ../data/square-disc.mesh -alpha 4.5 -o 3
                 mpirun -np 2 ex33p.py -m ../data/star.mesh -alpha 1.4 -o 3
                 mpirun -np 2 ex33p.py -m ../data/star.mesh -alpha 0.99 -o 3
                 mpirun -np 2 ex33p.py -m ../data/inline-quad.mesh -alpha 0.5 -o 3
                 mpirun -np 2 ex33p.py -m ../data/amr-quad.mesh -alpha 1.5 -o 3
                 mpirun -np 2 ex33p.py -m ../data/disc-nurbs.mesh -alpha 0.33 -o 3
                 mpirun -np 2 ex33p.py -m ../data/disc-nurbs.mesh -alpha 2.4 -o 3 -r 4
                 mpirun -np 2 ex33p.py -m ../data/l-shape.mesh -alpha 0.33 -o 3 -r 4
                 mpirun -np 2 ex33p.py -m ../data/l-shape.mesh -alpha 1.7 -o 3 -r 5
    Verification runs:
                 mpirun -np 2 ex33p.py -m ../data/inline-segment.mesh -ver -alpha 1.7 -o 2 -r 2
                 mpirun -np 2 ex33p.py -m ../data/inline-quad.mesh -ver -alpha 1.2 -o 2 -r 2
                 mpirun -np 2 ex33p.py -m ../data/amr-quad.mesh -ver -alpha 2.6 -o 2 -r 2
                 mpirun -np 2 ex33p.py -m ../data/inline-hex.mesh -ver -alpha 0.3 -o 2 -r 1
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


def run(alpha=0.5,
        order=1,
        num_refs=3,
        meshfile=None,
        visualization=False,
        verification=False,
        vef=False):

    coeffs = doubleArray()
    poles = doubleArray()
    progress_steps = 1

    power_of_laplace = floor(alpha)
    exponent_to_approximate = alpha - power_of_laplace

    if exponent_to_approximate > 1e-12:
        print("Approximating the fractional exponent " +
              str(exponent_to_approximate))

        from ex33_common import ComputePartialFractionApproximation
        poles, coeffs = ComputePartialFractionApproximation(
            exponent_to_approximate)

        # If the example is build without LAPACK, the exponent_to_approximate
        # might be modified by the function call above.
        alpha = exponent_to_approximate + power_of_laplace
        integer_order = False
    else:
        print("Treating integer order PDE.")
        integer_order = True

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    # 4. Refine the mesh to increase the resolution.
    for i in range(num_refs):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh

    # 5. Define a finite element space on the mesh.
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)

    fe_size = fespace.GlobalTrueVSize()
    if myid == 0:
        print("Number of finite element unknowns: " + str(fe_size))

    # 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    ess_tdof_list = intArray()
    if pmesh.bdr_attributes.Size() > 0:
        ess_bdr = intArray([1]*pmesh.bdr_attributes.Max())
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 7. Define diffusion coefficient, load, and solution GridFunction.
    @mfem.jit.scalar
    def f(x):
        val = 1.0
        for i in range(len(x)):
            val *= sin(pi*x[i])
        return (len(x)*pi**2)**alpha * val

    one = mfem.ConstantCoefficient(1.0)
    u = mfem.ParGridFunction(fespace)
    x = mfem.ParGridFunction(fespace)
    g = mfem.ParGridFunction(fespace)

    u.Assign(0.0)
    x.Assign(0.0)
    g.Assign(0.0)

    # 9. Set up the linear form b(.) for integer-order PDE solves.
    b = mfem.ParLinearForm(fespace)
    if verification:
        # This statement is only relevant for the verification of the code. It
        # uses a different f such that an analytic solution is known and easy
        # to compare with the numerical one. The FPDE becomes:
        # (-Δ)^α u = (2\pi ^2)^α sin(\pi x) sin(\pi y) on [0,1]^2
        # -> u(x,y) = sin(\pi x) sin(\pi y)
        b.AddDomainIntegrator(mfem.DomainLFIntegrator(f))
    else:
        b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))

    b.Assemble()

    # ------------------------------------------------------------------------
    # 10. Solve the PDE (-Δ)^N g = f, i.e. compute g = (-Δ)^{-1}^N f.
    # ------------------------------------------------------------------------

    if power_of_laplace > 0:
        # 10.1 Compute Stiffnes Matrix
        k = mfem.ParBilinearForm(fespace)
        k.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
        k.Assemble()

        # 10.2 Compute Mass Matrix
        m = mfem.ParBilinearForm(fespace)
        m.AddDomainIntegrator(mfem.MassIntegrator(one))
        m.Assemble()

        mass = mfem.HypreParMatrix()
        empty = intArray()
        m.FormSystemMatrix(empty, mass)

        #  10.3 Form the system of equations
        B = mfem.Vector()
        X = mfem.Vector()
        Op = mfem.OperatorPtr()
        k.FormLinearSystem(ess_tdof_list, g, b, Op, X, B)

        prec = mfem.HypreBoomerAMG()
        prec.SetPrintLevel(-1)
        cg = mfem.CGSolver(MPI.COMM_WORLD)
        cg.SetRelTol(1e-12)
        cg.SetMaxIter(2000)
        cg.SetPrintLevel(3)
        cg.SetPreconditioner(prec)
        cg.SetOperator(Op.Ptr())

        if myid == 0:
            print("\n".join(["",
                             "Computing (-Δ) ^ -" + str(power_of_laplace) + " ( f ) "]))

        for i in range(power_of_laplace):
            # 10.4 Solve the linear system Op X = B (N times).
            cg.Mult(B, X)

            # 10.5 Visualize the solution g of -Δ ^ N g = f in the last step
            if i == power_of_laplace - 1:
                # Needed for visualization and solution verification.
                k.RecoverFEMSolution(X, b, g)
                if integer_order and verification:
                    # For an integer order PDE, g is also our solution u.
                    u += g
                if visualization:
                    fout = mfem.socketstream("localhost", 19916)
                    fout.precision(8)
                    fout << "parallel " << num_procs << " " << myid << "\n"
                    fout << "solution\n" << pmesh << g

                    progress_steps = progress_steps + 1
                    title = ("Step " + str(progress_steps) + ": Solution of PDE -Δ ^ " +
                             str(power_of_laplace) + " g = f")
                    fout << "window_title '" + title + "'"
                    fout.flush()

            # 10.6 Prepare for next iteration (primal / dual space)
            mass.Mult(X, B)
            X.SetSubVectorComplement(ess_tdof_list, 0.0)

        # 10.7 Extract solution for the next step. The b now corresponds to the
        #      function g in the PDE.
        rm = mfem.fespace.GetRestrictionMatrix()
        rm.MultTranspose(B, b)

    # ------------------------------------------------------------------------
    # 11. Solve the fractional PDE by solving M integer order PDEs and adding
    #     up the solutions.
    # ------------------------------------------------------------------------
    if not integer_order:
        # Setup visualization.
        if visualization:
            xout = mfem.socketstream("localhost", 19916)
            uout = mfem.socketstream("localhost", 19916)

        # Iterate over all expansion coefficient that contribute to the
        # solution.
        for i in range(len(coeffs)):
            if myid == 0:
                print("\nSolving PDE -Δ u + " + "{:g}".format(-poles[i]) +
                      " u = " + "{:g}".format(coeffs[i]) + " g ")

            # 11.1 Reset GridFunction for integer-order PDE solve.
            x.Assign(0.0)

            # 11.2 Set up the bilinear form a(.,.) for integer-order PDE solve.
            a = mfem.ParBilinearForm(fespace)
            a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
            d_i = mfem.ConstantCoefficient(-poles[i])
            a.AddDomainIntegrator(mfem.MassIntegrator(d_i))
            a.Assemble()

            # 11.3 Assemble the bilinear form and the corresponding linear system.
            Op = mfem.OperatorPtr()
            B = mfem.Vector()
            X = mfem.Vector()
            a.FormLinearSystem(ess_tdof_list, x, b, Op, X, B)

            # 11.4 Solve the linear system A X = B.
            prec = mfem.HypreBoomerAMG()
            prec.SetPrintLevel(-1)
            cg = mfem.CGSolver(MPI.COMM_WORLD)
            cg.SetRelTol(1e-12)
            cg.SetMaxIter(2000)
            cg.SetPrintLevel(3)
            cg.SetPreconditioner(prec)
            cg.SetOperator(Op.Ptr())
            cg.Mult(B, X)

            # 11.5 Recover the solution as a finite element grid function.
            a.RecoverFEMSolution(X, b, x)

            # 11.6 Accumulate integer-order PDE solutions.
            x *= coeffs[i]
            u += x

            # 11.7 Send fractional PDE solution to a GLVis server.
            if visualization:
                oss_x = ("Step " + str(progress_steps) +
                         ": Solution of PDE -Δ u + " + "{:g}".format(-poles[i]) +
                         " u = " + "{:g}".format(coeffs[i]) + " g")
                xout << "parallel " << num_procs << " " << myid << "\n"
                xout << "solution\n" << pmesh << x
                xout << "window_title '" << oss_x << "'"
                xout.flush()

                oss_x = ("Step " + str(progress_steps+1) +
                         ": Solution of fractional PDE (-Δ)^ " + "{:g}".format(alpha) +
                         " u = f")
                uout << "parallel " << num_procs << " " << myid << "\n"
                uout << "solution\n" << pmesh << u
                uout << "window_title '" << oss_x << "'"
                uout.flush()

    # ------------------------------------------------------------------------
    # 12. (optional) Verify the solution.
    # ------------------------------------------------------------------------
    if verification:
        @mfem.jit.scalar
        def sol(x):
            val = 1.0
            for i in range(len(x)):
                val *= sin(pi*x[i])
            return val

        l2_error = u.ComputeL2Error(sol)

        if myid == 0:
            if dim == 1:
                analytic_solution = "sin(π x)"
                expected_mesh = "inline_segment.mesh"
            elif dim == 2:
                analytic_solution = "sin(π x) sin(π y)"
                expected_mesh = "inline_quad.mesh"
            else:
                analytic_solution = "sin(π x) sin(π y) sin(π z)"
                expected_mesh = "inline_hex.mesh"

            print("\n" + "="*80 +
                  "\n\nSolution Verification in " + str(dim) + "D \n\n" +
                  "Analytic solution : " + analytic_solution + "\n" +
                  "Expected mesh     : " + expected_mesh + "\n" +
                  "Your mesh         : " + meshfile + "\n" +
                  "L2 error          : " + "{:g}".format(l2_error) + "\n\n" +
                  "="*80)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex33 (Fractional PDE)')
    parser.add_argument('-m', '--mesh',
                        default="star.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space")
    parser.add_argument('-r', '--refs',
                        action='store', default=3, type=int,
                        help="Number of times to refine the mesh uniformly in serial")
    parser.add_argument("-alpha", "--alpha",
                        action='store',
                        type=float,
                        default=0.5,
                        help="Fractional exponent")
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        help='Disable GLVis visualization')
    parser.add_argument('-ver', '--verification',
                        action='store_true',
                        help="Use sinusoidal function (f) for analytic comparison.")

    try:
        from numba import jit
        HAS_NUMBA = True
    except ImportError:
        assert False, "This example requires numba to run"

    args = parser.parse_args()
    parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    visualization = False if args.no_visualization else True

    run(alpha=args.alpha,
        order=args.order,
        num_refs=args.refs,
        meshfile=meshfile,
        verification=args.verification,
        visualization=visualization,
        vef=args.verification)
