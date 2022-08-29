// Compile with: make meshOptMWE
// Sample run with gslib: meshOptMWE -m square01.mesh -o 2 -rs 2 -mid 80 -tid 5 -ni 50 -qo 4 -nor -vl 2 -ae 1
// Set -ae 0 for running without gslib


#include "mfem.hpp"
#include "common/mfem-common.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;


double discrete_size_2d(const Vector &x)
{
   int opt = 2;
   const double small = 0.001, big = 0.01;
   double val = 0.;

   const double xc = x(0) - 0.0, yc = x(1) - 0.5;
   const double r = sqrt(xc*xc + yc*yc);
   double r1 = 0.45; double r2 = 0.55; double sf=30.0;
   val = 0.5*(1+std::tanh(sf*(r-r1))) - 0.5*(1+std::tanh(sf*(r-r2)));

   val = std::max(0.,val);
   val = std::min(1.,val);

   return val * small + (1.0 - val) * big;
}

IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
IntegrationRules IntRulesCU(0, Quadrature1D::ClosedUniform);

int main(int argc, char *argv[])
{
   // 0. Set the method's default parameters.
   const char *mesh_file = "icf.mesh";
   int mesh_poly_deg     = 1;
   int rs_levels         = 0;
   int metric_id         = 1;
   int target_id         = 1;
   int quad_type         = 1;
   int quad_order        = 8;
   int solver_type       = 0;
   int solver_iter       = 20;
   double solver_rtol    = 1e-10;
   int solver_art_type   = 0;
   int lin_solver        = 2;
   int max_lin_iter      = 100;
   bool hradaptivity     = false;
   bool normalization    = false;
   bool visualization    = true;
   int verbosity_level   = 0;
   int adapt_eval        = 0;
   const char *devopt    = "cpu";
   int n_hr_iter         = 5;
   int n_h_iter          = 1;

   // 1. Parse command-line options.
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&mesh_poly_deg, "-o", "--order",
                  "Polynomial degree of mesh finite element space.");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&metric_id, "-mid", "--metric-id",
                  "Mesh optimization metric:\n\t"
                  "T-metrics\n\t"
                  "2  : 0.5|T|^2/tau-1                 -- 2D shape (condition number)\n\t"
                 );
   args.AddOption(&target_id, "-tid", "--target-id",
                  "Target (ideal element) type:\n\t"
                  "5: Ideal shape, given size (in physical space)");
   args.AddOption(&quad_type, "-qt", "--quad-type",
                  "Quadrature rule type:\n\t"
                  "1: Gauss-Lobatto\n\t"
                  "2: Gauss-Legendre\n\t"
                  "3: Closed uniform points");
   args.AddOption(&quad_order, "-qo", "--quad_order",
                  "Order of the quadrature rule.");
   args.AddOption(&solver_type, "-st", "--solver-type",
                  " Type of solver: (default) 0: Newton, 1: LBFGS");
   args.AddOption(&solver_iter, "-ni", "--newton-iters",
                  "Maximum number of Newton iterations.");
   args.AddOption(&solver_rtol, "-rtol", "--newton-rel-tolerance",
                  "Relative tolerance for the Newton solver.");
   args.AddOption(&solver_art_type, "-art", "--adaptive-rel-tol",
                  "Type of adaptive relative linear solver tolerance:\n\t"
                  "0: None (default)\n\t"
                  "1: Eisenstat-Walker type 1\n\t"
                  "2: Eisenstat-Walker type 2");
   args.AddOption(&lin_solver, "-ls", "--lin-solver",
                  "Linear solver:\n\t"
                  "0: l1-Jacobi\n\t"
                  "1: CG\n\t"
                  "2: MINRES\n\t"
                  "3: MINRES + Jacobi preconditioner\n\t"
                  "4: MINRES + l1-Jacobi preconditioner");
   args.AddOption(&max_lin_iter, "-li", "--lin-iter",
                  "Maximum number of iterations in the linear solve.");
   args.AddOption(&hradaptivity, "-hr", "--hr-adaptivity", "-no-hr",
                  "--no-hr-adaptivity",
                  "Enable hr-adaptivity.");
   args.AddOption(&normalization, "-nor", "--normalization", "-no-nor",
                  "--no-normalization",
                  "Make all terms in the optimization functional unitless.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&verbosity_level, "-vl", "--verbosity-level",
                  "Set the verbosity level - 0, 1, or 2.");
   args.AddOption(&adapt_eval, "-ae", "--adaptivity-evaluator",
                  "0 - Advection based (DEFAULT), 1 - GSLIB.");
   args.AddOption(&n_hr_iter, "-nhr", "--n_hr_iter",
                  "Number of hr-adaptivity iterations.");
   args.AddOption(&n_h_iter, "-nh", "--n_h_iter",
                  "Number of h-adaptivity iterations per r-adaptivity"
                  "iteration.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   if (hradaptivity)
   {
      MFEM_VERIFY(strcmp(devopt,"cpu")==0, "HR-adaptivity is currently only"
                  " supported on cpus.");
   }

   // 2. Initialize and refine the starting mesh.
   Mesh *mesh = new Mesh(mesh_file, 1, 1, false);
   for (int lev = 0; lev < rs_levels; lev++) { mesh->UniformRefinement(); }
   const int dim = mesh->Dimension();

   FiniteElementCollection *fec = new H1_FECollection(mesh_poly_deg, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim);

   mesh->SetNodalFESpace(fespace);

   Vector b(0);

   GridFunction x(fespace);
   mesh->SetNodalGridFunction(&x);
   x.SetTrueVector();
   x.SetFromTrueVector();

   // 9. Save the starting (prior to the optimization) mesh to a file. This
   //    output can be viewed later using GLVis: "glvis -m perturbed.mesh".
   {
      ofstream mesh_ofs("perturbed.mesh");
      mesh->Print(mesh_ofs);
   }

   // 10. Store the starting (prior to the optimization) positions.
   GridFunction x0(fespace);
   x0 = x;

   TMOP_QualityMetric *metric = new TMOP_Metric_080(0.5);

   TargetConstructor::TargetType target_t;
   TargetConstructor *target_c = NULL;
   H1_FECollection ind_fec(mesh_poly_deg, dim);
   FiniteElementSpace ind_fes(mesh, &ind_fec);
   GridFunction size(&ind_fes);

   switch (target_id)
   {
      case 5: // Discrete size 2D or 3D
      {
         target_t = TargetConstructor::IDEAL_SHAPE_GIVEN_SIZE;
         DiscreteAdaptTC *tc = new DiscreteAdaptTC(target_t);
         if (adapt_eval == 0)
         {
            tc->SetAdaptivityEvaluator(new AdvectorCG());
         }
         else
         {
#ifdef MFEM_USE_GSLIB
            tc->SetAdaptivityEvaluator(new InterpolatorFP);
#else
            MFEM_ABORT("MFEM is not built with GSLIB.");
#endif
         }
         if (dim == 2)
         {
            FunctionCoefficient size_coeff(discrete_size_2d);
            size.ProjectCoefficient(size_coeff);
         }
         else
         {
             MFEM_ABORT("only dim == 2 supported for this MWE.");
         }
         tc->SetSerialDiscreteTargetSize(size);
         target_c = tc;
         break;
      }
      default: cout << "Unknown target_id: " << target_id << endl; return 3;
   }
   if (target_c == NULL)
   {
      target_c = new TargetConstructor(target_t);
   }
   target_c->SetNodes(x0);
   TMOP_Integrator *tmop_integ = new TMOP_Integrator(metric, target_c);

   // Setup the quadrature rules for the TMOP integrator.
   IntegrationRules *irules = NULL;
   switch (quad_type)
   {
      case 1: irules = &IntRulesLo; break;
      case 2: irules = &IntRules; break;
      case 3: irules = &IntRulesCU; break;
      default: cout << "Unknown quad_type: " << quad_type << endl; return 3;
   }
   tmop_integ->SetIntegrationRules(*irules, quad_order);

   if (normalization) { tmop_integ->EnableNormalization(x0); }

   NonlinearForm a(fespace);
   ConstantCoefficient *metric_coeff1 = NULL;
   a.AddDomainIntegrator(tmop_integ);

   // For HR tests, the energy is normalized by the number of elements.
   const double init_energy = a.GetGridFunctionEnergy(x);

   // Visualize the starting mesh and metric values.
   // Note that for combinations of metrics, this only shows the first metric.
   if (visualization)
   {
      char title[] = "Initial metric values";
      vis_tmop_metric_s(mesh_poly_deg, *metric, *target_c, *mesh, title, 0);
   }

   // 13. Fix all boundary nodes, or fix only a given component depending on the
   //     boundary attributes of the given mesh. Attributes 1/2/3 correspond to
   //     fixed x/y/z components of the node. Attribute 4 corresponds to an
   //     entirely fixed node. Other boundary attributes do not affect the node
   //     movement boundary conditions.
   Array<int> vdofs;
   Array<int> ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 1;
   a.SetEssentialBC(ess_bdr);

   // 14. As we use the Newton method to solve the resulting nonlinear system,
   //     here we setup the linear solver for the system's Jacobian.
   Solver *S = NULL, *S_prec = NULL;
   const double linsol_rtol = 1e-12;
   if (lin_solver == 0)
   {
      S = new DSmoother(1, 1.0, max_lin_iter);
   }
   else if (lin_solver == 1)
   {
      CGSolver *cg = new CGSolver;
      cg->SetMaxIter(max_lin_iter);
      cg->SetRelTol(linsol_rtol);
      cg->SetAbsTol(0.0);
      cg->SetPrintLevel(verbosity_level >= 2 ? 3 : -1);
      S = cg;
   }
   else
   {
      MINRESSolver *minres = new MINRESSolver;
      minres->SetMaxIter(max_lin_iter);
      minres->SetRelTol(linsol_rtol);
      minres->SetAbsTol(0.0);
      if (verbosity_level > 2) { minres->SetPrintLevel(1); }
      minres->SetPrintLevel(verbosity_level == 2 ? 3 : -1);
      if (lin_solver == 3 || lin_solver == 4)
      {
        auto ds = new DSmoother((lin_solver == 3) ? 0 : 1, 1.0, 1);
        ds->SetPositiveDiagonal(true);
        S_prec = ds;
         minres->SetPreconditioner(*S_prec);
      }
      S = minres;
   }

   // Perform the nonlinear optimization.
   const IntegrationRule &ir =
      irules->Get(fespace->GetFE(0)->GetGeomType(), quad_order);
   TMOPNewtonSolver solver(ir, solver_type);
   solver.SetIntegrationRules(*irules, quad_order);
   if (solver_type == 0)
   {
      // Specify linear solver when we use a Newton-based solver.
      solver.SetPreconditioner(*S);
   }
   solver.SetMaxIter(solver_iter);
   solver.SetRelTol(solver_rtol);
   solver.SetAbsTol(0.0);
   if (solver_art_type > 0)
   {
      solver.SetAdaptiveLinRtol(solver_art_type, 0.5, 0.9);
   }
   solver.SetPrintLevel(verbosity_level >= 1 ? 1 : -1);

   TMOPHRSolver hr_solver(*mesh, a, solver,
                          x, false, hradaptivity,
                          mesh_poly_deg, metric_id,
                          n_hr_iter, n_h_iter);
   hr_solver.AddGridFunctionForUpdate(&x0);
   hr_solver.Mult();

   // 15. Save the optimized mesh to a file. This output can be viewed later
   //     using GLVis: "glvis -m optimized.mesh".
   {
      ofstream mesh_ofs("optimized.mesh");
      mesh_ofs.precision(14);
      mesh->Print(mesh_ofs);
   }

   const double fin_energy = a.GetGridFunctionEnergy(x);
   std::cout << std::scientific << std::setprecision(4);
   cout << "Initial strain energy: " << init_energy << endl;
   cout << "  Final strain energy: " << fin_energy << endl;
   cout << "The strain energy decreased by: "
        << (init_energy - fin_energy) * 100.0 / init_energy << " %." << endl;

   // 16. Visualize the final mesh and metric values.
   if (visualization)
   {
      char title[] = "Final metric values";
      vis_tmop_metric_s(mesh_poly_deg, *metric, *target_c, *mesh, title, 600);
   }

   // 17. Visualize the mesh displacement.
   if (visualization)
   {
      osockstream sock(19916, "localhost");
      sock << "solution\n";
      mesh->Print(sock);
      x0 -= x;
      x0.Save(sock);
      sock.send();
      sock << "window_title 'Displacements'\n"
           << "window_geometry "
           << 1200 << " " << 0 << " " << 600 << " " << 600 << "\n"
           << "keys jRmclA" << endl;
   }

   delete S;
   delete S_prec;
   delete metric_coeff1;
   delete target_c;
   delete metric;
   delete fespace;
   delete fec;
   delete mesh;

   return 0;
}
