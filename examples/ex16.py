'''
   MFEM example 16

   How to run:
      python <arguments>

   Example of arguments:
      ex16.py -m inline-tri.mesh
      ex16.py -m disc-nurbs.mesh -tf 2
      ex16.py -s 1 -a 0.0 -k 1.0
      ex16.py -s 2 -a 1.0 -k 0.0
      ex16.py -s 3 -a 0.5 -k 0.5 -o 4
      ex16.py -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
      ex16.py -m fichera-q2.mesh
      ex16.py -m escher.mesh
      ex16.py -m beam-tet.mesh -tf 10 -dt 0.1
      ex16.py -m amr-quad.mesh -o 4 -r 0
      ex16.py -m amr-hex.mesh -o 2 -r 0

'''
import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from mfem import path

import mfem.ser as mfem
from mfem.ser import intArray


class ConductionOperator(mfem.PyTimeDependentOperator):
    def __init__(self, fespace, alpha, kappa, u):
        mfem.PyTimeDependentOperator.__init__(
            self, fespace.GetTrueVSize(), 0.0)
        rel_tol = 1e-8
        self.alpha = alpha
        self.kappa = kappa
        self.T = None
        self.K = None
        self.M = None
        self.fespace = fespace

        self.ess_tdof_list = intArray()
        self.Mmat = mfem.SparseMatrix()
        self.Kmat = mfem.SparseMatrix()
        self.M_solver = mfem.CGSolver()
        self.M_prec = mfem.DSmoother()
        self.T_solver = mfem.CGSolver()
        self.T_prec = mfem.DSmoother()
        self.z = mfem.Vector(self.Height())

        self.M = mfem.BilinearForm(fespace)
        self.M.AddDomainIntegrator(mfem.MassIntegrator())
        self.M.Assemble()
        self.M.FormSystemMatrix(self.ess_tdof_list, self.Mmat)

        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(rel_tol)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(30)
        self.M_solver.SetPrintLevel(0)
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.Mmat)

        self.T_solver.iterative_mode = False
        self.T_solver.SetRelTol(rel_tol)
        self.T_solver.SetAbsTol(0.0)
        self.T_solver.SetMaxIter(100)
        self.T_solver.SetPrintLevel(0)
        self.T_solver.SetPreconditioner(self.T_prec)

        self.SetParameters(u)

    def Mult(self, u, u_dt):
        # Compute:
        #  du_dt = M^{-1}*-K(u) for du_dt
        self.Kmat.Mult(u, self.z)
        self.z.Neg()   # z = -z
        self.M_solver.Mult(self.z, du_dt)

    def ImplicitSolve(self, dt, u, du_dt):
        # Solve the equation:
        #    du_dt = M^{-1}*[-K(u + dt*du_dt)]
        #    for du_dt
        if self.T is None:
            self.T = mfem.Add(1.0, self.Mmat, dt, self.Kmat)
            current_dt = dt
            self.T_solver.SetOperator(self.T)
        self.Kmat.Mult(u, self.z)
        self.z.Neg()
        self.T_solver.Mult(self.z, du_dt)

    def SetParameters(self, u):
        u_alpha_gf = mfem.GridFunction(self.fespace)
        u_alpha_gf.SetFromTrueDofs(u)
        for i in range(u_alpha_gf.Size()):
            u_alpha_gf[i] = self.kappa + self.alpha * u_alpha_gf[i]

        self.K = mfem.BilinearForm(self.fespace)
        u_coeff = mfem.GridFunctionCoefficient(u_alpha_gf)
        self.K.AddDomainIntegrator(mfem.DiffusionIntegrator(u_coeff))
        self.K.Assemble()
        self.K.FormSystemMatrix(self.ess_tdof_list, self.Kmat)
        self.T = None


class InitialTemperature(mfem.PyCoefficient):
    def EvalValue(self, x):
        xx = np.array(x)
        norm2 = np.sqrt(np.sum(xx**2))
        if norm2 < 0.5:
            return 2.0
        return 1.0


parser = ArgParser(description='Ex16')
parser.add_argument('-m', '--mesh',
                    default='star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-r', '--refine',
                    action='store', default=2, type=int,
                    help="Number of times to refine the mesh uniformly, -1 for auto.")
parser.add_argument('-o', '--order',
                    action='store', default=2, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-s', '--ode-solver',
                    action='store', default=3, type=int,
                    help='\n'.join(["ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3",
                                    "\t\t 11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4."]))
parser.add_argument('-t', '--t-final',
                    action='store', default=0.5, type=float,
                    help="Final time; start time is 0.")
parser.add_argument("-dt", "--time-step",
                    action='store', default=0.01, type=float,
                    help="Time step.")
parser.add_argument('-a', '--alpha',
                    action='store', default=0.01, type=float,
                    help='Alpha coefficient')
parser.add_argument('-k', '--kappa',
                    action='store', default=0.5, type=float,
                    help='Kappa coefficient')
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')
parser.add_argument('-vs', '--visualization-steps',
                    action='store', default=5, type=int,
                    help="Visualize every n-th timestep.")

args = parser.parse_args()
ref_levels = args.refine
order = args.order
dt = args.time_step
t_final = args.t_final
alpha = args.alpha
kappa = args.kappa
visualization = args.visualization
vis_steps = args.visualization_steps
ode_solver_type = args.ode_solver
meshfile = expanduser(
    join(dirname(__file__), '..', 'data', args.mesh))
parser.print_options(args)

# 2. Read the mesh from the given mesh file. We can handle triangular,
#    quadrilateral, tetrahedral and hexahedral meshes with the same code.
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()

# 3. Define the ODE solver used for time integration. Several implicit
#    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
#    explicit Runge-Kutta methods are available.
if ode_solver_type == 1:
    ode_solver = BackwardEulerSolver()
elif ode_solver_type == 2:
    ode_solver = mfem.SDIRK23Solver(2)
elif ode_solver_type == 3:
    ode_solver = mfem.SDIRK33Solver()
elif ode_solver_type == 11:
    ode_solver = ForwardEulerSolver()
elif ode_solver_type == 12:
    ode_solver = mfem.RK2Solver(0.5)
elif ode_solver_type == 13:
    ode_solver = mfem.RK3SSPSolver()
elif ode_solver_type == 14:
    ode_solver = mfem.RK4Solver()
elif ode_solver_type == 22:
    ode_solver = mfem.ImplicitMidpointSolver()
elif ode_solver_type == 23:
    ode_solver = mfem.SDIRK23Solver()
elif ode_solver_type == 24:
    ode_solver = mfem.SDIRK34Solver()
else:
    print("Unknown ODE solver type: " + str(ode_solver_type))
    exit
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

# 6. Set the initial conditions for u. All boundaries are considered
#    natural.
u_0 = InitialTemperature()
u_gf.ProjectCoefficient(u_0)
u = mfem.Vector()
u_gf.GetTrueDofs(u)

# 7. Initialize the conduction operator and the visualization.
oper = ConductionOperator(fespace, alpha, kappa, u)
u_gf.SetFromTrueDofs(u)

mesh.Print('ex16.mesh', 8)
u_gf.Save('ex16-init.gf', 8)

if visualization:
    sout = mfem.socketstream("localhost", 19916)
    sout.precision(8)
    sout << "solution\n" << mesh << u_gf
    sout << "pause\n"
    sout.flush()
    print("GLVis visualization paused.")
    print(" Press space (in the GLVis window) to resume it.")

# 8. Perform time-integration (looping over the time iterations, ti, with a
#    time-step dt).
ode_solver.Init(oper)
t = 0.0
ti = 1
last_step = False
while not last_step:
    if (t + dt >= t_final - dt/2):
        last_step = True

    t, dt = ode_solver.Step(u, t, dt)

    if (last_step or (ti % vis_steps) == 0):
        # if True:
        print("step " + str(ti) + ", t = " + "{:g}".format(t))
        u_gf.SetFromTrueDofs(u)
        if (visualization):
            sout << "solution\n" << mesh << u_gf
            sout.flush()
    ti = ti + 1
    oper.SetParameters(u)

# 9. Save the final solution. This output can be viewed later using GLVis:
#    "glvis -m ex16.mesh -g ex16-final.gf".
u_gf.Save('ex16-final.gf', 8)
