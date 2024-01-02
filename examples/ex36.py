'''
   PyMFEM example 36

   See c++ version in the MFEM library for more detail

   Sample runs:  python ex36.py -o 2
                 python ex36.py -o 2 -r 4


'''
from numpy import sqrt

visualization = True

def run(refs=args.refs,
        order=args.order,
        max_it=args.max_it,
        step=args.step,
        alpha=args.alpha)


    # 2. Read the mesh from the mesh file.
    mesh_file = "../data/disc-nurbs.mesh";    
    mesh_file = expanduser(
        join(os.path.dirname(__file__), mesh_file.split("/"))
    mesh = mfem.Mesh(mesh_file, 1, 1)       
    dim = mesh.Dimension();

    # 3. Postprocess the mesh.
    # 3A. Refine the mesh to increase the resolution.
    for l in range(refs):
        mesh.UniformRefinement();

    # 3B. Interpolate the geometry after refinement to control geometry error.
    # NOTE: Minimum second-order interpolation is used to improve the accuracy.
    curvature_order = max((order,2))
    mesh.SetCurvature(curvature_order)

    # 3C. Rescale the domain to a unit circle (radius = 1).
    nodes = mesh.GetNodes()
    scale = 2*sqrt(2)
    nodes = node/scale

    # 4. Define the necessary finite element spaces on the mesh.
    H1fec = mfem.H1_FECollection(order+1, dim)
    H1fes = mfem.FiniteElementSpace(mesh, H1fec)

    L2fec = mfem.L2_FECollection(order-1, dim)
    L2fes = mfem.FiniteElementSpace(mesh, L2fec)

    print("Number of H1 finite element unknowns: " +
          str(H1fes.GetTrueVSize()))
    print("Number of L2 finite element unknowns: " + 
          str(L2fes.GetTrueVSize()))

    offset = mfem.intArray((0, H1fes.GetVSize(), L2fes.GetVSize()))
    offsets.PartialSum()

    x = mfem.BlockVector(offsets)
    rhs = mfem.BlockVector(offsets)
    x.Assign(0.0)
    rhs.Assign(0.0)

    # 5. Determine the list of true (i.e. conforming) essential boundary dofs.
    if mesh.bdr_attributes.Size() > 0:
        ess_bdr = mfem.intArray([1]*mesh.bdr_attributes.Max())
    else:
        ess_bdr = mfem.intArray()

    # 6. Define an initial guess for the solution.
    def IC_func(x):
        r0 = 1.0;
        rr = 0.0;
        for i in range(x.Size()):
            rr += x[i]*x[i]
        return r0*r0 - rr

    one = mfem.ConstantCoefficient(1.0)
    zero = mfem.ConstantCoefficient(0.0)        

   # 7. Define the solution vectors as a finite element grid functions
   #    corresponding to the fespaces.
   u_gf = mfem.GridFunction()
   delta_psi_gf = mfem.GridFunction()
   u_gf.MakeRef(H1fes,x,offsets[0])
   delta_psi_gf.MakeRef&L2fes,x,offsets[1])
   delta_psi_gf.Assign(0.0)

   u_old_gf = mfem.GridFunction(H1fes)
   psi_old_gf = mfem.GridFunction(L2fes)
   psi_gf = mfem.GridFunction(L2fes)
   u_old_gf.Assign(0.0)
   psi_old_gf.Assign(0.0)

   # 8. Define the function coefficients for the solution and use them to
   #    initialize the initial guess
   exact_coef = mfem.FunctionCoefficient(exact_solution_obstacle)
   exact_grad_coef = mfem.VectorFunctionCoefficient(dim,
                                            exact_solution_gradient_obstacle)
   IC_coef = mfem.FunctionCoefficient(IC_func)
   f = mfem.ConstantCoefficient(0.0)
   obstacle = mfem.FunctionCoefficient(spherical_obstacle)
   u_gf.ProjectCoefficient(IC_coef)
   u_old_gf.GetDataArray()[:] = u_gf.GetDataArray()

   # 9. Initialize the slack variable ψₕ = ln(uₕ)
   ln_u = mfem.LogarithmGridFunctionCoefficient(u_gf, obstacle)
   psi_gf.ProjectCoefficient(ln_u)
   psi_old_gf.GetDataArray()[:] = psi_gf.GetDataArray()

   if visualization:
       sol_sock = mfem.socketstream("localhost", 19916)
       sol_sock.precision(8)

   # 10. Iterate
   total_iterations = 0
   increment_u = 0.
   for k in range(max_it):
      u_tmp = mfem.GridFunction(&1fes)
      u_tmp.GetDataArray()[:] = u_old_gf.GetDataArray()

      print("\nOUTER ITERATION " + str(k+1))

      for 0 in range(10):
         total_iterations += 1

         alpha_cf = mfem.ConstantCoefficient(alpha)

         b0 = mfem.LinearForm()
         b1 = mfem.LinearForm()         
         b0.Update(H1fes, rhs.GetBlock(0), 0)
         b1.Update(L2fes, rhs.GetBlock(1), 0)

         exp_psi = mfem.ExponentialGridFunctionCoefficient (psi_gf, zero);
         ProductCoefficient neg_exp_psi(-1.0,exp_psi);
         GradientGridFunctionCoefficient grad_u_old(&u_old_gf);
         ProductCoefficient alpha_f(alpha, f);
         GridFunctionCoefficient psi_cf(&psi_gf);
         GridFunctionCoefficient psi_old_cf(&psi_old_gf);
         SumCoefficient psi_old_minus_psi(psi_old_cf, psi_cf, 1.0, -1.0);

         b0.AddDomainIntegrator(mfem.DomainLFIntegrator(alpha_f))
         b0.AddDomainIntegrator(mfem.DomainLFIntegrator(psi_old_minus_psi))
         b0.Assemble();

         b1.AddDomainIntegrator(mfem.DomainLFIntegrator(exp_psi))
         b1.AddDomainIntegrator(mfem.DomainLFIntegrator(obstacle))
         b1.Assemble();

         BilinearForm a00(&H1fes);
         a00.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
         a00.AddDomainIntegrator(new DiffusionIntegrator(alpha_cf));
         a00.Assemble();
         a00.EliminateEssentialBC(ess_bdr,x.GetBlock(0),rhs.GetBlock(0),
                                  mfem::Operator::DIAG_ONE);
         a00.Finalize();
         SparseMatrix &A00 = a00.SpMat();

         MixedBilinearForm a10(&H1fes,&L2fes);
         a10.AddDomainIntegrator(new MixedScalarMassIntegrator());
         a10.Assemble();
         a10.EliminateTrialDofs(ess_bdr, x.GetBlock(0), rhs.GetBlock(1));
         a10.Finalize();
         SparseMatrix &A10 = a10.SpMat();

         SparseMatrix *A01 = Transpose(A10);

         BilinearForm a11(&L2fes);
         a11.AddDomainIntegrator(new MassIntegrator(neg_exp_psi));
         // NOTE: Shift the spectrum of the Hessian matrix for additional
         //       stability (Quasi-Newton).
         ConstantCoefficient eps_cf(-1e-6);
         if (order == 1)
         {
            // NOTE: ∇ₕuₕ = 0 for constant functions.
            //       Therefore, we use the mass matrix to shift the spectrum
            a11.AddDomainIntegrator(new MassIntegrator(eps_cf));
         }
         else
         {
            a11.AddDomainIntegrator(new DiffusionIntegrator(eps_cf));
         }
         a11.Assemble();
         a11.Finalize();
         SparseMatrix &A11 = a11.SpMat();

         BlockOperator A(offsets);
         A.SetBlock(0,0,&A00);
         A.SetBlock(1,0,&A10);
         A.SetBlock(0,1,A01);
         A.SetBlock(1,1,&A11);

         BlockDiagonalPreconditioner prec(offsets);
         prec.SetDiagonalBlock(0,new GSSmoother(A00));
         prec.SetDiagonalBlock(1,new GSSmoother(A11));
         prec.owns_blocks = 1;

         GMRES(A,prec,rhs,x,0,10000,500,1e-12,0.0);

         u_gf.MakeRef(&H1fes, x.GetBlock(0), 0);
         delta_psi_gf.MakeRef(&L2fes, x.GetBlock(1), 0);

         u_tmp -= u_gf;
         double Newton_update_size = u_tmp.ComputeL2Error(zero);
         u_tmp = u_gf;

         double gamma = 1.0;
         delta_psi_gf *= gamma;
         psi_gf += delta_psi_gf;

         if (visualization)
         {
            sol_sock << "solution\n" << mesh << u_gf << "window_title 'Discrete solution'"
                     << flush;
            mfem::out << "Newton_update_size = " << Newton_update_size << endl;
         }

         delete A01;

         if (Newton_update_size < increment_u)
         {
            break;
         }
      }

      u_tmp = u_gf;
      u_tmp -= u_old_gf;
      increment_u = u_tmp.ComputeL2Error(zero);

      mfem::out << "Number of Newton iterations = " << j+1 << endl;
      mfem::out << "Increment (|| uₕ - uₕ_prvs||) = " << increment_u << endl;

      u_old_gf = u_gf;
      psi_old_gf = psi_gf;

      if (increment_u < tol || k == max_it-1)
      {
         break;
      }

      double H1_error = u_gf.ComputeH1Error(&exact_coef,&exact_grad_coef);
      mfem::out << "H1-error  (|| u - uₕᵏ||)       = " << H1_error << endl;

   }

   mfem::out << "\n Outer iterations: " << k+1
             << "\n Total iterations: " << total_iterations
             << "\n Total dofs:       " << H1fes.GetTrueVSize() + L2fes.GetTrueVSize()
             << endl;

   // 11. Exact solution.
   if (visualization)
   {
      socketstream err_sock(vishost, visport);
      err_sock.precision(8);

      GridFunction error_gf(&H1fes);
      error_gf.ProjectCoefficient(exact_coef);
      error_gf -= u_gf;

      err_sock << "solution\n" << mesh << error_gf << "window_title 'Error'"  <<
               flush;
   }

   {
      double L2_error = u_gf.ComputeL2Error(exact_coef);
      double H1_error = u_gf.ComputeH1Error(&exact_coef,&exact_grad_coef);

      ExponentialGridFunctionCoefficient u_alt_cf(psi_gf,obstacle);
      GridFunction u_alt_gf(&L2fes);
      u_alt_gf.ProjectCoefficient(u_alt_cf);
      double L2_error_alt = u_alt_gf.ComputeL2Error(exact_coef);

      mfem::out << "\n Final L2-error (|| u - uₕ||)          = " << L2_error <<
                endl;
      mfem::out << " Final H1-error (|| u - uₕ||)          = " << H1_error << endl;
      mfem::out << " Final L2-error (|| u - ϕ - exp(ψₕ)||) = " << L2_error_alt <<
                endl;
   }

   return 0;
}

double LogarithmGridFunctionCoefficient::Eval(ElementTransformation &T,
                                              const IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "grid function is not set");

   double val = u->GetValue(T, ip) - obstacle->Eval(T, ip);
   return max(min_val, log(val));
}

double ExponentialGridFunctionCoefficient::Eval(ElementTransformation &T,
                                                const IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "grid function is not set");

   double val = u->GetValue(T, ip);
   return min(max_val, max(min_val, exp(val) + obstacle->Eval(T, ip)));
}

double spherical_obstacle(const Vector &pt)
{
   double x = pt(0), y = pt(1);
   double r = sqrt(x*x + y*y);
   double r0 = 0.5;
   double beta = 0.9;

   double b = r0*beta;
   double tmp = sqrt(r0*r0 - b*b);
   double B = tmp + b*b/tmp;
   double C = -b/tmp;

   if (r > b)
   {
      return B + r * C;
   }
   else
   {
      return sqrt(r0*r0 - r*r);
   }
}

double exact_solution_obstacle(const Vector &pt)
{
   double x = pt(0), y = pt(1);
   double r = sqrt(x*x + y*y);
   double r0 = 0.5;
   double a =  0.348982574111686;
   double A = -0.340129705945858;

   if (r > a)
   {
      return A * log(r);
   }
   else
   {
      return sqrt(r0*r0-r*r);
   }
}

void exact_solution_gradient_obstacle(const Vector &pt, Vector &grad)
{
   double x = pt(0), y = pt(1);
   double r = sqrt(x*x + y*y);
   double r0 = 0.5;
   double a =  0.348982574111686;
   double A = -0.340129705945858;

   if (r > a)
   {
      grad(0) =  A * x / (r*r);
      grad(1) =  A * y / (r*r);
   }
   else
   {
      grad(0) = - x / sqrt( r0*r0 - r*r );
      grad(1) = - y / sqrt( r0*r0 - r*r );
   }
}




if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex36(Energy Minimization)')
   

    parser.add_argument("-o", "--order",
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).")
    parser.add_argument('-r', '--refs', 
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument('-mi', '--max-it', 
                        action='store', default=1, type=int,
                        help="Maximum number of iterations")
    help = "\n".join(("Stopping criteria based on the difference between",
                      "successive solution updates"))
    parser.add_argument('-tol', '--tol',
                        action='store', default=1e-5, type=float,                        
                        help=help)
    parser.add_argument('-step', '--step',
                        action='store', default=1.0, type=float,                        
                        help="Step size alpha")
    parser.add_argument("-no-vis","--no-visualization",
                        action='store_false', default=True,
                        help='Enable GLVis visualization')

    
    args = parser.parse_args()
    parser.print_options(args)

    globals["visualization"] = args.no_visualization

    run(refs=args.refs,
        order=args.order,
        max_it=args.max_it,
        step=args.step,
        alpha=args.alpha)

