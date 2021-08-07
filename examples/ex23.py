'''
   MFEM example 23
      See c++ version in the MFEM library for more detail 
'''
import os
import io
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array, floor, log, sum


def run(mesh_file="",
        ref_levels=2,
        order=2,
        ode_solver_type=10,
        t_final=0.5,
        dt=1e-2,
        speed=1.0,
        dirichlet=True,
        visit=True,
        visualization=True,
        vis_steps=5):

    class WaveOperator(mfem.SecondOrderTimeDependentOperator):
        def __init__(self, fespace, ess_bdr, speed):
            mfem.SecondOrderTimeDependentOperator.__init__(
                self, fespace.GetTrueVSize(), 0.0)

            self.ess_tdof_list = mfem.intArray()

            rel_tol = 1e-8
            fespace.GetEssentialTrueDofs(ess_bdr, self.ess_tdof_list)

            c2 = mfem.ConstantCoefficient(speed*speed)
            K = mfem.BilinearForm(fespace)
            K.AddDomainIntegrator(mfem.DiffusionIntegrator(c2))
            K.Assemble()

            self.Kmat0 = mfem.SparseMatrix()
            self.Kmat = mfem.SparseMatrix()
            dummy = mfem.intArray()
            K.FormSystemMatrix(dummy, self.Kmat0)
            K.FormSystemMatrix(self.ess_tdof_list, self.Kmat)
            self.K = K

            self.Mmat = mfem.SparseMatrix()
            M = mfem.BilinearForm(fespace)
            M.AddDomainIntegrator(mfem.MassIntegrator())
            M.Assemble()
            M.FormSystemMatrix(self.ess_tdof_list, self.Mmat)
            self.M = M

            M_solver = mfem.CGSolver()
            M_prec = mfem.DSmoother()
            M_solver.iterative_mode = False
            M_solver.SetRelTol(rel_tol)
            M_solver.SetAbsTol(0.0)
            M_solver.SetMaxIter(30)
            M_solver.SetPrintLevel(0)
            M_solver.SetPreconditioner(M_prec)
            M_solver.SetOperator(self.Mmat)
            self.M_prec = M_prec
            self.M_solver = M_solver

            T_solver = mfem.CGSolver()
            T_prec = mfem.DSmoother()
            T_solver.iterative_mode = False
            T_solver.SetRelTol(rel_tol)
            T_solver.SetAbsTol(0.0)
            T_solver.SetMaxIter(100)
            T_solver.SetPrintLevel(0)
            T_solver.SetPreconditioner(T_prec)
            self.T_prec = T_prec
            self.T_solver = T_solver
            self.T = None

        def Mult(self, u, du_dt, d2udt2):
            # Compute:
            #    d2udt2 = M^{-1}*-K(u)
            # for d2udt2
            z = mfem.Vector(u.Size())
            self.Kmat.Mult(u, z)
            z.Neg()  # z = -z
            self.M_solver.Mult(z, d2udt2)

        def ImplicitSolve(self, fac0, fac1, u, dudt, d2udt2):
            # Solve the equation:
            #    d2udt2 = M^{-1}*[-K(u + fac0*d2udt2)]
            # for d2udt2
            if self.T is None:
                self.T = mfem.Add(1.0, self.Mmat, fac0, self.Kmat)
                self.T_solver.SetOperator(self.T)
            z = mfem.Vector(u.Size())
            self.Kmat0.Mult(u, z)
            z.Neg()

            # iterate over Array<int> :D
            for j in self.ess_tdof_list:
                z[j] = 0.0

            self.T_solver.Mult(z, d2udt2)

        def SetParameters(self, u):
            self.T = None

    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()

    # 3. Define the ODE solver used for time integration. Several second order
    #    time integrators are available.
    if ode_solver_type <= 10:
        ode_solver = mfem.GeneralizedAlpha2Solver(ode_solver_type/10)
    elif ode_solver_type == 11:
        ode_solver = mfem.AverageAccelerationSolver()
    elif ode_solver_type == 12:
        ode_solver = mfem.LinearAccelerationSolver()
    elif ode_solver_type == 13:
        ode_solver = mfem.CentralDifferenceSolver()
    elif ode_solver_type == 14:
        ode_solver = mfem.FoxGoodwinSolver()
    else:
        print("Unknown ODE solver type: " + str(ode_solver_type))

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement, where 'ref_levels' is a
    #    command-line parameter.
    for lev in range(ref_levels):
        mesh.UniformRefinement()

    # 5. Define the vector finite element space representing the current and the
    #    initial temperature, u_ref.
    fe_coll = mfem.H1_FECollection(order, dim)
    fespace = mfem.FiniteElementSpace(mesh, fe_coll)

    fe_size = fespace.GetTrueVSize()
    print("Number of temperature unknowns: " + str(fe_size))

    u_gf = mfem.GridFunction(fespace)
    dudt_gf = mfem.GridFunction(fespace)

    # 6. Set the initial conditions for u. All boundaries are considered
    #    natural.
    class cInitialSolution(mfem.PyCoefficient):
        def EvalValue(self, x):
            norm2 = sum(x**2)
            return exp(-norm2*30)

    class cInitialRate(mfem.PyCoefficient):
        def EvalValue(self, x):
            return 0

    u_0 = cInitialSolution()
    dudt_0 = cInitialRate()

    u_gf.ProjectCoefficient(u_0)
    u = mfem.Vector()
    u_gf.GetTrueDofs(u)

    dudt_gf.ProjectCoefficient(dudt_0)
    dudt = mfem.Vector()
    dudt_gf.GetTrueDofs(dudt)

    # 7. Initialize the conduction operator and the visualization.
    ess_bdr = mfem.intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr.SetSize(mesh.bdr_attributes.Max())
        if (dirichlet):
            ess_bdr.Assign(1)
        else:
            ess_bdr.Assigne(0)

    oper = WaveOperator(fespace, ess_bdr, speed)

    u_gf.SetFromTrueDofs(u)

    mesh.Print("ex23.mesh", 8)
    output = io.StringIO()
    output.precision = 8
    u_gf.Save(output)
    dudt_gf.Save(output)
    fid = open("ex23-init.gf", 'w')
    fid.write(output.getvalue())
    fid.close()

    if visit:
        visit_dc = mfem.VisItDataCollection("Example23", mesh)
        visit_dc.RegisterField("solution", u_gf)
        visit_dc.RegisterField("rate", dudt_gf)
        visit_dc.SetCycle(0)
        visit_dc.SetTime(0.0)
        visit_dc.Save()

    if visualization:
        sout = mfem.socketstream("localhost", 19916)
        if not sout.good():
            print("Unable to connect to GLVis server at localhost:19916")
            visualization = False
            print("GLVis visualization disabled.")
        else:
            sout.precision(precision)
            sout << "solution\n" << mesh << dudt_gf
            sout << "pause\n"
            sout.flush()
            print(
                "GLVis visualization paused. Press space (in the GLVis window) to resume it.")

    # 8. Perform time-integration (looping over the time iterations, ti, with a
    #    time-step dt).
    ode_solver.Init(oper)
    t = 0.0

    last_step = False
    ti = 0
    while not last_step:
        ti += 1
        if t + dt >= t_final - dt/2:
            last_step = True

        t, dt = ode_solver.Step(u, dudt, t, dt)

        if last_step or (ti % vis_steps == 0):
            print("step " + str(ti) + ", t = " + "{:g}".format(t))

            u_gf.SetFromTrueDofs(u)
            dudt_gf.SetFromTrueDofs(dudt)
            if visualization:
                sout << "solution\n" << mesh << u_gf
                sout.flush()

            if visit:
                visit_dc.SetCycle(ti)
                visit_dc.SetTime(t)
                visit_dc.Save()

        oper.SetParameters(u)

    # 9. Save the final solution. This output can be viewed later using GLVis:
    #    "glvis -m ex23.mesh -g ex23-final.gf".
    output = io.StringIO()
    output.precision = 8
    u_gf.Save(output)
    dudt_gf.Save(output)
    fid = open("ex23-final.gf", "w")
    fid.write(output.getvalue())
    fid.close()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser
    parser = ArgParser(description="Ex23 (Wave problem)")
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-r', '--refine',
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly before parallel")
    parser.add_argument('-o', '--order',
                        action='store', default=2, type=int,
                        help="Finite element order (polynomial degree)")
    help_ode = '\n'.join(["ODE solver: [0--10] \t- GeneralizedAlpha(0.1 * s),",
                          "11 \t - Average Acceleration,",
                          "12 \t - Linear Acceleration",
                          "13 \t- CentralDifference",
                          "14 \t- FoxGoodwin"])
    parser.add_argument('-s', '--ode-solver',
                        action='store', default=10, type=int,
                        help=help_ode)
    parser.add_argument('-tf', '--t-final',
                        action='store', default=0.5, type=float,
                        help="Final time; start time is 0.")
    parser.add_argument('-dt', '--time-step',
                        action='store', default=1e-2, type=float,
                        help="Time step")
    parser.add_argument("-c", "--speed",
                        action='store', default=1.0, type=float,
                        help="Wave speed.")
    parser.add_argument("-neu",  "--neumann",
                        action='store_true', default=False,
                        help="BC switch.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true', default=True,
                        help='Enable GLVis visualization')
    parser.add_argument('-visit', '--visit-datafiles',
                        action='store_true', default=True,
                        help="Save data files for VisIt (visit.llnl.gov) visualization.")
    parser.add_argument("-vs", "--visualization-steps",
                        action='store', default=5,  type=int,
                        help="Visualize every n-th timestep.")

    args = parser.parse_args()
    parser.print_options(args)
    mesh_file = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    run(mesh_file=mesh_file,
        ref_levels=args.refine,
        order=args.order,
        ode_solver_type=args.ode_solver,
        t_final=args.t_final,
        dt=args.time_step,
        speed=args.speed,
        dirichlet=(not args.neumann),
        visit=args.visit_datafiles,
        vis_steps=args.visualization_steps,
        visualization=args.visualization)
