#
#   based on meshOptMWE.cpp from Ketan Mittal
#
# python test_tmop.py -m ../data/square01.mesh -o 2 -rs 2 -mid 80 -tid 5 -ni 50 -qo 4 -vl 2 -ae 0
import sys
import os
from os.path import expanduser, join
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == '-p':
    import mfem.par as mfem
    use_parallel = True
    from mfem.common.mpi_debug import nicePrint as print
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank
    sys.argv = [sys.argv[0]] + sys.argv[2:]

else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0

class discrete_size_2d(mfem.PyCoefficient):
    def EvalValue(self, x):
        opt = 2;
        small = 0.001
        big = 0.01
        val = 0.;
        
        xc = x[0] - 0.0
        yc = x[1] - 0.5
        r = np.sqrt(xc*xc + yc*yc)
        r1 = 0.45
        r2 = 0.55
        sf=30.0
        val = 0.5*(1+np.tanh(sf*(r-r1))) - 0.5*(1+np.tanh(sf*(r-r2)))

        val = max(0.,val)
        val = min(1.,val)

        return val * small + (1.0 - val) * big
    
def run(args):
    mesh_file = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    mesh_poly_deg = args.order
    rs_levels = args.refine_serial
    metric_id = args.metric_id
    target_id = args.target_id
    quad_type = args.quad_type
    quad_order = args.quad_order
    solver_type = args.solver_type
    lin_solver = args.lin_solver
    normalization = args.normalization

    verbosity_level = args.verbosity_level
    
    n_hr_iter=args.n_hr_iter
    n_h_iter=args.n_h_iter
   
    solver_iter = args.newton_iters
    solver_rtol= args.newton_rel_tolerance
    solver_art_type = args.adaptive_rel_tol
    max_lin_iter = args.lin_iter
    hradaptivity = args.hr_adaptivity
    visualization = not args.no_visualization
    adapt_eval= args.adaptivity_evaluator

    devopt    = "cpu";

    # 2. Initialize and refine the starting mesh.
    mesh = mfem.Mesh(mesh_file, 1, 1, False)
    for i in range(rs_levels):
        mesh.UniformRefinement()
    dim = mesh.Dimension();

    fec = mfem.H1_FECollection(mesh_poly_deg, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec, dim)

    mesh.SetNodalFESpace(fespace)

    b = mfem.Vector(0)

    x = mfem.GridFunction(fespace)
    mesh.SetNodalGridFunction(x)
    x.SetTrueVector()
    x.SetFromTrueVector()
                              
    # 9. Save the starting (prior to the optimization) mesh to a file. This
    #    output can be viewed later using GLVis: "glvis -m perturbed.mesh".
    mesh.Print("perturbed.mesh")

    # 10. Store the starting (prior to the optimization) positions.
    x0 = mfem.GridFunction(fespace)
    x0.Assign(x)

    metric = mfem.tmop.TMOP_Metric_080(0.5)

    ind_fec = mfem.H1_FECollection(mesh_poly_deg, dim)
    ind_fes = mfem.FiniteElementSpace(mesh, ind_fec)
    size = mfem.GridFunction(ind_fes)

    
    if target_id == 5: # Discrete size 2D or 3D
         target_t = mfem.tmop.TargetConstructor.IDEAL_SHAPE_GIVEN_SIZE
                              
         tc = mfem.tmop.DiscreteAdaptTC(target_t)
                              
         if adapt_eval == 0:
             tc.SetAdaptivityEvaluator(mfem.tmop.AdvectorCG())
         else:
             if "InterpolatorFP" in dir(mfem.tmop):
                 evaluator = mfem.tmop.InterpolatorFP()
                 tc.SetAdaptivityEvaluator(evaluator)
             else:
                 assert False, "MFEM is not built with GSLIB."
         if dim == 2:
             #size_coeff = mfem.FunctionCoefficient(discrete_size_2d)
             size_coeff = discrete_size_2d()
             size.ProjectCoefficient(size_coeff)
         else:
             assert False, "only dim == 2 supported for this MWE."
                              
         tc.SetSerialDiscreteTargetSize(size)
         target_c = tc;

    else:
        print("Unknown target_id: " + str(target_id))
        return

    if target_c is None:
       target_c = mfem.tmop.TargetConstructor(target_t);

    target_c.SetNodes(x0)
    
    tmop_integ = mfem.tmop.TMOP_Integrator(metric, target_c)
 
    # Setup the quadrature rules for the TMOP integrator.
    if quad_type == 1:
        irules = mfem.IntegrationRules(0, mfem.Quadrature1D.GaussLobatto)
    elif quad_type == 2:        
        irules = mfem.IntRules
    elif quad_type == 3:        
        irules = mfem.IntegrationRules(0, mfem.Quadrature1D.ClosedUniform)
    else:
        print( "Unknown quad_type: " + str(quad_type))
        return 3

    tmop_integ.SetIntegrationRules(irules, quad_order)

    if normalization:
        tmop_integ.EnableNormalization(x0)

    a = mfem.NonlinearForm(fespace)
    a.AddDomainIntegrator(tmop_integ)

    # For HR tests, the energy is normalized by the number of elements.
    init_energy = a.GetGridFunctionEnergy(x);
    
    # Visualize the starting mesh and metric values.
    # Note that for combinations of metrics, this only shows the first metric.
    if visualization:
        mfem.tmop.vis_tmop_metric_s(mesh_poly_deg, metric, target_c, mesh, "Initial metric values", 0);

    # 13. Fix all boundary nodes, or fix only a given component depending on the
    #     boundary attributes of the given mesh. Attributes 1/2/3 correspond to
    #     fixed x/y/z components of the node. Attribute 4 corresponds to an
    #     entirely fixed node. Other boundary attributes do not affect the node
    #     movement boundary conditions.
    ess_bdr = mfem.intArray([1]*mesh.bdr_attributes.Max())
    a.SetEssentialBC(ess_bdr);


    # 14. As we use the Newton method to solve the resulting nonlinear system,
    #     here we setup the linear solver for the system's Jacobian.
    linsol_rtol = 1e-12;
    if lin_solver == 0:
        S = mfem.DSmoother(1, 1.0, max_lin_iter)
    elif lin_solver == 1:
        cg = mfem.CGSolver()
        cg.SetMaxIter(max_lin_iter)
        cg.SetRelTol(linsol_rtol)
        cg.SetAbsTol(0.0)
        cg.SetPrintLevel(3 if verbosity_level >= 2 else -1)
        S = cg
    else:
        minres = mfem.MINRESSolver()
        minres.SetMaxIter(max_lin_iter)
        minres.SetRelTol(linsol_rtol)
        minres.SetAbsTol(0.0)
        if verbosity_level > 2:
            minres.SetPrintLevel(1)
        minres.SetPrintLevel(3 if verbosity_level == 2 else -1)
        if lin_solver == 3 or lin_solver == 4:
             ds = mfem.DSmoother((0 if lin_solver == 3 else 1), 1.0, 1)
             ds.SetPositiveDiagonal(True)
             minres.SetPreconditioner(ds)
        S = minres;
   
    #/ Perform the nonlinear optimization.
    ir = irules.Get(fespace.GetFE(0).GetGeomType(), quad_order)
    solver = mfem.tmop.TMOPNewtonSolver(ir, solver_type)
    solver.SetIntegrationRules(irules, quad_order)
    if solver_type == 0:
        # Specify linear solver when we use a Newton-based solver.
        solver.SetPreconditioner(S)

    print(dir(solver))
    solver.SetMaxIter(solver_iter)
    solver.SetRelTol(solver_rtol)
    solver.SetAbsTol(0.0)
    if solver_art_type > 0:
        solver.SetAdaptiveLinRtol(solver_art_type, 0.5, 0.9)

    solver.SetPrintLevel(1 if verbosity_level >= 1 else -1)

    hr_solver = mfem.tmop.TMOPHRSolver(mesh, a, solver, x, False, hradaptivity,
                                  mesh_poly_deg, metric_id, n_hr_iter, n_h_iter)
    hr_solver.AddGridFunctionForUpdate(x0)
    hr_solver.Mult()

    # 15. Save the optimized mesh to a file. This output can be viewed later
    #     using GLVis: "glvis -m optimized.mesh".
    mesh.Print("optimized.mesh", 14)

    fin_energy = a.GetGridFunctionEnergy(x)
    print("Initial strain energy: " + "{:g}".format(init_energy))
    print("  Final strain energy: " + "{:g}".format(fin_energy))
    print("The strain energy decreased by: " + 
         "{:g}".format((init_energy - fin_energy) * 100.0 / init_energy))

    # 16. Visualize the final mesh and metric values.
    if visualization:
       mfem.tmop.vis_tmop_metric_s(mesh_poly_deg, metric, target_c, mesh, "Final metric values", 600);

    # 17. Visualize the mesh displacement.
    if visualization:
        x0 -= x
        sock = mfem.socketstream("localhost", 19916)
        sock << "solution\n" << mesh << x0
        sock.flush()
        sock << "window_title 'Displacements'\n" << "window_geometry "
        sock << 1200 << " " << 0 << " " << 600 << " " << 600 << "\n"
        sock << "keys jRmclA"
        sock.flush()             


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='meshOptMWE')
    parser.add_argument('-m', '--mesh',
                        default='square01.mesh',  # icf.mesh
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")
    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly in serial")
    parser.add_argument("-mid", "--metric-id",
                        action="store", default=80, type=int,
                        help="\n".join(["Mesh optimization metric:",
                                        "\tT-metrics",
                                        "2  : 0.5|T|^2/tau-1                 -- 2D shape (condition number)\n\t"]))
    parser.add_argument("-tid", "--target-id",
                        action="store", default=5, type=int,                        
                        help="\n".join(["Target (ideal element) type:",
                                        "\t5: Ideal shape, given size (in physical space)"]))
    parser.add_argument("-qt", "--quad-type",
                        action="store", default=1, type=int,                                                
                        help="\n".join(["Quadrature rule type:",
                                        "\t1: Gauss-Lobatto",
                                        "\t2: Gauss-Legendre"
                                        "\t3: Closed uniform points"]))
    parser.add_argument("-qo", "--quad_order",
                        action="store", default=4, type=int,                                                
                         help="Order of the quadrature rule.")
    parser.add_argument("-st", "--solver-type",
                        action="store", default=0, type=int,
                        help = " Type of solver: (default) 0: Newton, 1: LBFGS")
    parser.add_argument("-ni", "--newton-iters",
                        action="store", default=80, type=int,                        
                        help="Maximum number of Newton iterations.")
    parser.add_argument("-rtol", "--newton-rel-tolerance",
                        action="store", default=1e-10, type=float,
                        help="Relative tolerance for the Newton solver.")
    parser.add_argument("-art", "--adaptive-rel-tol",
                        action="store", default=0, type=int,                        
                        help="\n".join(["Type of adaptive relative linear solver tolerance:",
                                        "\t0: None (default)",
                                        "\t1: Eisenstat-Walker type 1",
                                        "\t2: Eisenstat-Walker type 2"]))
    parser.add_argument("-ls", "--lin-solver",
                        action="store", default=2, type=int,                                                 
                        help="\n".join(["Linear solver:",
                                        "\t0: l1-Jacobi",
                                        "\t1: CG",
                                        "\t2: MINRES",
                                        "\t3: MINRES + Jacobi preconditioner",
                                        "\t4: MINRES + l1-Jacobi preconditioner"]))
    parser.add_argument("-li", "--lin-iter",
                        action="store", default=100, type=int,     
                        help="Maximum number of iterations in the linear solve.")
    parser.add_argument("-hr", "--hr-adaptivity", 
                        action='store_true',
                        help="Enable hr-adaptivity.")
    parser.add_argument("-nor", "--normalization", 
                        action='store_true',
                        help="Make all terms in the optimization functional unitless.")
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument("-vl", "--verbosity-level",
                        action="store", default=2, type=int,                                                 
                        help="Set the verbosity level - 0, 1, or 2.")


    parser.add_argument("-ae", "--adaptivity-evaluator",
                        action="store", default=1, type=int,                                                 
                        help="0 - Advection based (DEFAULT), 1 - GSLIB.");

    parser.add_argument("-nhr", "--n_hr_iter",
                        action="store", default=5, type=int,
                        help="Number of hr-adaptivity iterations.")

    parser.add_argument("-nh", "--n_h_iter",
                        action="store", default=1, type=int,
                        help="Number of h-adaptivity iterations per r-adaptivity")
                  
    args = parser.parse_args()
    parser.print_options(args)
                  
    run(args)
