'''
   MFEM example 16

   How to run:
      mpirun -np 4 python <arguments>

   Example of arguments:
      ex16p.py -m inline-tri.mesh
      ex16p.py -m disc-nurbs.mesh -tf 2
      ex16p.py -s 1 -a 0.0 -k 1.0
      ex16p.py -s 2 -a 1.0 -k 0.0
      ex16p.py -s 3 -a 0.5 -k 0.5 -o 4
      ex16p.py -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
      ex16p.py -m fichera-q2.mesh
      ex16p.py -m escher.mesh
      ex16p.py -m beam-tet.mesh -tf 10 -dt 0.1
      ex16p.py -m amr-quad.mesh -o 4 -rs 0 -rp 0
      ex16p.py -m amr-hex.mesh -o 2 -rs 0 -rp 0
'''
import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

from mfem.par import intArray, add_vector
import mfem.par as mfem
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


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
        self.Mmat = mfem.HypreParMatrix()
        self.Kmat = mfem.HypreParMatrix()
        self.M_solver = mfem.CGSolver(fespace.GetComm())
        self.M_prec = mfem.HypreSmoother()
        self.T_solver = mfem.CGSolver(fespace.GetComm())
        self.T_prec = mfem.HypreSmoother()
        self.z = mfem.Vector(self.Height())

        self.M = mfem.ParBilinearForm(fespace)
        self.M.AddDomainIntegrator(mfem.MassIntegrator())
        self.M.Assemble()
        self.M.FormSystemMatrix(self.ess_tdof_list, self.Mmat)

        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(rel_tol)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(100)
        self.M_solver.SetPrintLevel(0)
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)
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
        u_alpha_gf = mfem.ParGridFunction(self.fespace)
        u_alpha_gf.SetFromTrueDofs(u)
        for i in range(u_alpha_gf.Size()):
            u_alpha_gf[i] = self.kappa + self.alpha * u_alpha_gf[i]

        self.K = mfem.ParBilinearForm(self.fespace)
        u_coeff = mfem.GridFunctionCoefficient(u_alpha_gf)
        self.K.AddDomainIntegrator(mfem.DiffusionIntegrator(u_coeff))
        self.K.Assemble(0)
        self.K.FormSystemMatrix(self.ess_tdof_list, self.Kmat)
        self.T = None


class InitialTemperature(mfem.PyCoefficient):
    def EvalValue(self, x):
        xx = np.array(x)
        norm2 = np.sqrt(float(np.sum(xx**2)))
        if norm2 < 0.5:
            return 2.0
        return 1.0


parser = ArgParser(description='Ex16p')
parser.add_argument('-m', '--mesh',
                    default='star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-rs', '--refine-serial',
                    action='store', default=2, type=int,
                    help="Number of times to refine the mesh uniformly in serial")
parser.add_argument('-rp', '--refine-parallel',
                    action='store', default=1, type=int,
                    help="Number of times to refine the mesh uniformly in parallel")

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
ser_ref_levels = args.refine_serial
par_ref_levels = args.refine_parallel
order = args.order
dt = args.time_step
t_final = args.t_final
alpha = args.alpha
kappa = args.kappa
visualization = args.visualization
vis_steps = args.visualization_steps
ode_solver_type = args.ode_solver
if myid == 0:
    parser.print_options(args)

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

# 3. Read the serial mesh from the given mesh file on all processors. We can
#    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
#    with the same code.
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()

# 4. Define the ODE solver used for time integration. Several implicit
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
# 5. Refine the mesh in serial to increase the resolution. In this example
#    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
#    a command-line parameter.
for lev in range(ser_ref_levels):
    mesh.UniformRefinement()

# 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
#    this mesh further in parallel to increase the resolution. Once the
#    parallel mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
for x in range(par_ref_levels):
    pmesh.UniformRefinement()

# 7. Define the vector finite element space representing the current and the
#    initial temperature, u_ref.
fe_coll = mfem.H1_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fe_coll)

fe_size = fespace.GlobalTrueVSize()
if myid == 0:
    print("Number of temperature unknowns: " + str(fe_size))
u_gf = mfem.ParGridFunction(fespace)

# 8. Set the initial conditions for u. All boundaries are considered
#    natural.
u_0 = InitialTemperature()
u_gf.ProjectCoefficient(u_0)
u = mfem.Vector()
u_gf.GetTrueDofs(u)

# 9. Initialize the conduction operator and the visualization.
oper = ConductionOperator(fespace, alpha, kappa, u)
u_gf.SetFromTrueDofs(u)

smyid = '{:0>6d}'.format(myid)
mesh_name = "ex16-mesh."+smyid
sol_name = "ex16-init."+smyid
pmesh.Print(mesh_name, 8)
u_gf.Save(sol_name, 8)

if visualization:
    sout = mfem.socketstream("localhost", 19916)
    sout.send_text("parallel " + str(num_procs) + " " + str(myid))
    sout.precision(8)
    sout.send_solution(pmesh, u_gf)
    sout.send_text("pause")
    sout.flush()
    if myid == 0:
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
        if myid == 0:
            print("step " + str(ti) + ", t = " + "{:g}".format(t))
        u_gf.SetFromTrueDofs(u)
        if (visualization):
            sout.send_text("parallel " + str(num_procs) + " " + str(myid))
            sout.send_solution(pmesh, u_gf)
            sout.flush()
    ti = ti + 1
    oper.SetParameters(u)

# 11. Save the final solution in parallel. This output can be viewed later
#     using GLVis: "glvis -np <np> -m ex16-mesh -g ex16-final".
sol_name = "ex16-final."+smyid
u_gf.Save(sol_name, 8)
