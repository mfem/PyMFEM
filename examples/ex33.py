'''
   MFEM example 33

   See c++ version in the MFEM library for more detail

   Sample runs:  python ex33.py -m ../data/square-disc.mesh -alpha 0.33 -o 2
                 python ex33.py -m ../data/square-disc.mesh -alpha 4.5 -o 3
                 python ex33.py -m ../data/star.mesh -alpha 1.4 -o 3
                 python ex33.py -m ../data/star.mesh -alpha 0.99 -o 3
                 python ex33.py -m ../data/inline-quad.mesh -alpha 0.5 -o 3
                 python ex33.py -m ../data/amr-quad.mesh -alpha 1.5 -o 3
                 python ex33.py -m ../data/disc-nurbs.mesh -alpha 0.33 -o 3
                 python ex33.py -m ../data/disc-nurbs.mesh -alpha 2.4 -o 3 -r 4
                 python ex33.py -m ../data/l-shape.mesh -alpha 0.33 -o 3 -r 4
                 python ex33.py -m ../data/l-shape.mesh -alpha 1.7 -o 3 -r 5
    Verification runs:
                 python ex33.py -m ../data/inline-segment.mesh -ver -alpha 1.7 -o 2 -r 2
                 python ex33.py -m ../data/inline-quad.mesh -ver -alpha 1.2 -o 2 -r 2
                 python ex33.py -m ../data/amr-quad.mesh -ver -alpha 2.6 -o 2 -r 2
                 python ex33.py -m ../data/inline-hex.mesh -ver -alpha 0.3 -o 2 -r 1
'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor


def run(alpha=0.5,
        order=1,
        rs=3,
        meshfile=None,
        visualization=False,
        vef=False):

   coeffs = doubleArray()
   poles = doubleArray()
   progress_steps = 1;

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
       alpha = exponent_to_approximate + power_of_laplace;
       integer_order = False
   else:
       print("Treating integer order PDE.")
       integer_order = True;

   print("coeffs", coeffs)
   print("poles", poles)

   mesh = mfem.Mesh(mesh_file, 1, 1)
   dim = mesh.Dimension()

   # 4. Refine the mesh to increase the resolution.
   for i in range(num_refs):
      mesh.UniformRefinement()

   # 5. Define a finite element space on the mesh.
   fec = mfem.H1_FECollection(order, dim)
   fespace = mfem.FiniteElementSpace(mesh, fec)

   print("Number of finite element unknowns: " + str(fespace.GetTrueVSize())

   # 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   ess_tdof_list=intArray()
   if mesh.bdr_attributes.Size() > 0:
       ess_bdr=intAttay([1]*mesh.bdr_attributes.Max())
       fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

   # 7. Define diffusion coefficient, load, and solution GridFunction.
   @ mfem.jit.scalar
   def f(x):
       val=1.0
       for i in range(x.Size()):
           val *= sin(pi*x[i])
       return (len(x)*pi**2)**alpha * val

   one=mfem.ConstantCoefficient(1.0)
   u=mfem.GridFunction(fespace)
   x=mfem.GridFunction(fespace)
   g=mfem.GridFunction(fespace)

   u.Assign(0.0)
   x.Assign(0.0)
   g.Assign(0.0)

   # 9. Set up the linear form b(.) for integer-order PDE solves.
   b=mfem.LinearForm(fespace)
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

   // 8. Prepare for visualization.
   char vishost[]="localhost";
   int  visport=19916;

   # ------------------------------------------------------------------------
   # 10. Solve the PDE (-Δ)^N g = f, i.e. compute g = (-Δ)^{-1}^N f.
   # ------------------------------------------------------------------------

   if power_of_laplace > 0:
       # 10.1 Compute Stiffnes Matrix
       k=mfem.BilinearForm(fespace)
       k.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
       k.Assemble()

       # 10.2 Compute Mass Matrix
       m=mfem.BilinearForm(fespace)
       m.AddDomainIntegrator(mfem.MassIntegrator(one))
       m.Assemble()

       mass=mfem.SparseMatrix()
       empty=intArray()
       m.FormSystemMatrix(empty, mass)

       #  10.3 Form the system of equations
       B=mfem.Vector()
       X=mfem.Vector()
       Op=mfem.OperatorPtr()
       k.FormLinearSystem(ess_tdof_list, g, b, Op, X, B)

       M = mfem.GSSmoother(SparseMatrix &)(*Op));

       print("\n".join(["",
                        "Computing (-Δ) ^ -" + str(power_of_laplace) + " ( f ) "]))

       for i in range(power_of_laplace):
           # 10.4 Solve the linear system Op X = B (N times).
           mfem.PCG(Op, M, B, X, 3, 300, 1e-12, 0.0)

           # 10.5 Visualize the solution g of -Δ ^ N g = f in the last step
           if i == power_of_laplace - 1:
               # Needed for visualization and solution verification.
               k.RecoverFEMSolution(X, b, g)
               if integer_order and verification:
                    # For an integer order PDE, g is also our solution u.
                    u += g
               if (visualization)
                   fout = mfem.socketstream("localhost", 19916)
                   fout.precision(8)
                   fout << "solution\n" << mesh << g

                   progress_steps = progress_steps + 1
                   title = ("Step " + str(progress_steps) + ": Solution of PDE -Δ ^ " +
                            str(power_of_laplace) + " g = f")
                   fout << "window_title '" + title "'"
                   fout.flush()

           # 10.6 Prepare for next iteration (primal / dual space)
           mass.Mult(X, B)
           X.SetSubVectorComplement(ess_tdof_list, 0.0)

       # 10.7 Extract solution for the next step. The b now corresponds to the
       #      function g in the PDE.
       R = mfem.fespace.GetRestrictionMatrix()
       if R is not None:
           R.MultTranspose(B, b)
       else:
           b.Assign(B)

   # ------------------------------------------------------------------------
   # 11. Solve the fractional PDE by solving M integer order PDEs and adding
   #     up the solutions.
   # ------------------------------------------------------------------------
   if not integer_order:
      # Setup visualization.
      if visualization:
          xout = mfem.socketstream("localhost", 19916)
          uout = mfem.socketstream("localhost", 19916)
         xout.open(vishost, visport);
         xout.precision(8);
         uout.open(vishost, visport);
         uout.precision(8);
      }
      // Iterate over all expansion coefficient that contribute to the
      // solution.
      for (int i = 0; i < coeffs.Size(); i++)
      {
         mfem::out << "\nSolving PDE -Δ u + " << -poles[i]
                   << " u = " << coeffs[i] << " g " << endl;


         // 11.1 Reset GridFunction for integer-order PDE solve.
         x = 0.0;

         // 11.2 Set up the bilinear form a(.,.) for integer-order PDE solve.
         BilinearForm a(&fespace);
         a.AddDomainIntegrator(new DiffusionIntegrator(one));
         ConstantCoefficient d_i(-poles[i]);
         a.AddDomainIntegrator(new MassIntegrator(d_i));
         a.Assemble();

         // 11.3 Assemble the bilinear form and the corresponding linear system.
         OperatorPtr A;
         Vector B, X;
         a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

         // 11.4 Solve the linear system A X = B.
         GSSmoother M((SparseMatrix&)(*A));

         PCG(*A, M, B, X, 3, 300, 1e-12, 0.0);

         // 11.5 Recover the solution as a finite element grid function.
         a.RecoverFEMSolution(X, b, x);

         // 11.6 Accumulate integer-order PDE solutions.
         x *= coeffs[i];
         u += x;

         // 11.7 Send fractional PDE solution to a GLVis server.
         if (visualization)
         {
            oss_x.str(""); oss_x.clear();
            oss_x << "Step " << progress_steps
                  << ": Solution of PDE -Δ u + " << -poles[i]
                  << " u = " << coeffs[i] << " g";
            xout << "solution\n" << mesh << x
                 << "window_title '" << oss_x.str() << "'" << flush;

            oss_u.str(""); oss_u.clear();
            oss_u << "Step " << progress_steps + 1
                  << ": Solution of fractional PDE (-Δ)^" << alpha
                  << " u = f";
            uout << "solution\n" << mesh << u
                 << "window_title '" << oss_u.str() << "'"
                 << flush;
         }
      }
   }

   // ------------------------------------------------------------------------
   // 12. (optional) Verify the solution.
   // ------------------------------------------------------------------------
   if (verification)
   {
      auto solution = [] (const Vector &x)
      {
         double val = 1.0;
         for (int i=0; i<x.Size(); i++)
         {
            val *= sin(M_PI*x(i));
         }
         return val;
      };
      FunctionCoefficient sol(solution);
      double l2_error = u.ComputeL2Error(sol);

      string analytic_solution,expected_mesh;
      switch (dim)
      {
         case 1:
            analytic_solution = "sin(π x)";
            expected_mesh = "inline_segment.mesh";
            break;
         case 2:
            analytic_solution = "sin(π x) sin(π y)";
            expected_mesh = "inline_quad.mesh";
            break;
         default:
            analytic_solution = "sin(π x) sin(π y) sin(π z)";
            expected_mesh = "inline_hex.mesh";
            break;
      }

      mfem::out << "\n" << string(80,'=')
                << "\n\nSolution Verification in "<< dim << "D \n\n"
                << "Analytic solution : " << analytic_solution << "\n"
                << "Expected mesh     : " << expected_mesh <<"\n"
                << "Your mesh         : " << mesh_file << "\n"
                << "L2 error          : " << l2_error << "\n\n"
                << string(80,'=') << endl;
   }

   return 0;
}




    
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
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
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

    order = args.order
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    alpha = args.alpha


    run(alpha=alpha,
        order=order,
        rs=args.refs,
        meshfile=meshfile,
        visualization=visualization,
        vef=args.verification)

    

    
