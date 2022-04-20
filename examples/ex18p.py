'''
   MFEM example 18 parallel
      This is a version of Example 18 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
import mfem.par as mfem

from ex18_common import FE_Evolution, InitialCondition, RiemannSolver, DomainIntegrator, FaceIntegrator
from mfem.common.arg_parser import ArgParser

from os.path import expanduser, join, dirname
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

# Equation constant parameters.(using globals to share them with ex18_common)
import ex18_common

# 1. Initialize MPI.from mpi4py import MPI

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

parser = ArgParser(description='Ex18p')
parser.add_argument('-m', '--mesh',
                    default='periodic-square.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-p', '--problem',
                    action='store', default=1, type=int,
                    help='Problem setup to use. See options in velocity_function().')
parser.add_argument('-rs', '--refine_serial',
                    action='store', default=0, type=int,
                    help="Number of times to refine the mesh uniformly before parallel.")
parser.add_argument('-rp', '--refine_parallel',
                    action='store', default=1, type=int,
                    help="Number of times to refine the mesh uniformly after parallel.")
parser.add_argument('-o', '--order',
                    action='store', default=3, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-s', '--ode_solver',
                    action='store', default=4, type=int,
                    help="ODE solver: 1 - Forward Euler,\n\t" +
                    "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6.")
parser.add_argument('-tf', '--t_final',
                    action='store', default=2.0, type=float,
                    help="Final time; start time is 0.")
parser.add_argument("-dt", "--time_step",
                    action='store', default=-0.01, type=float,
                    help="Time step.")
parser.add_argument('-c', '--cfl_number',
                    action='store', default=0.3, type=float,
                    help="CFL number for timestep calculation.")
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')
parser.add_argument('-vs', '--visualization-steps',
                    action='store', default=50, type=float,
                    help="Visualize every n-th timestep.")

args = parser.parse_args()
mesh = args.mesh
ser_ref_levels = args.refine_serial
par_ref_levels = args.refine_parallel
order = args.order
ode_solver_type = args.ode_solver
t_final = args.t_final
dt = args.time_step
cfl = args.cfl_number
visualization = args.visualization
vis_steps = args.visualization_steps

if myid == 0:
    parser.print_options(args)

device = mfem.Device('cpu')
if myid == 0:            
    device.Print()

ex18_common.num_equation = 4
ex18_common.specific_heat_ratio = 1.4
ex18_common.gas_constant = 1.0
ex18_common.problem = args.problem
num_equation = ex18_common.num_equation


# 3. Read the mesh from the given mesh file. This example requires a 2D
#    periodic mesh, such as ../data/periodic-square.mesh.
meshfile = expanduser(join(dirname(__file__), '..', 'data', mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()

# 4. Define the ODE solver used for time integration. Several explicit
#    Runge-Kutta methods are available.
ode_solver = None
if ode_solver_type == 1:
    ode_solver = mfem.ForwardEulerSolver()
elif ode_solver_type == 2:
    ode_solver = mfem.RK2Solver(1.0)
elif ode_solver_type == 3:
    ode_solver = mfem.RK3SSolver()
elif ode_solver_type == 4:
    ode_solver = mfem.RK4Solver()
elif ode_solver_type == 6:
    ode_solver = mfem.RK6Solver()
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
for lev in range(par_ref_levels):
    pmesh.UniformRefinement()

# 7. Define the discontinuous DG finite element space of the given
#    polynomial order on the refined mesh.
fec = mfem.DG_FECollection(order, dim)
# Finite element space for a scalar (thermodynamic quantity)
fes = mfem.ParFiniteElementSpace(pmesh, fec)
# Finite element space for a mesh-dim vector quantity (momentum)
dfes = mfem.ParFiniteElementSpace(pmesh, fec, dim, mfem.Ordering.byNODES)
# Finite element space for all variables together (total thermodynamic state)
vfes = mfem.ParFiniteElementSpace(
    pmesh, fec, num_equation, mfem.Ordering.byNODES)

assert fes.GetOrdering() == mfem.Ordering.byNODES, "Ordering must be byNODES"
glob_size = vfes.GlobalTrueVSize()
if myid == 0:
    print("Number of unknowns: " + str(glob_size))

# 8. Define the initial conditions, save the corresponding mesh and grid
#    functions to a file. This can be opened with GLVis with the -gc option.
#    The solution u has components {density, x-momentum, y-momentum, energy}.
#    These are stored contiguously in the BlockVector u_block.

offsets = [k*vfes.GetNDofs() for k in range(num_equation+1)]
offsets = mfem.intArray(offsets)
u_block = mfem.BlockVector(offsets)

#  Momentum grid function on dfes for visualization.
mom = mfem.ParGridFunction(dfes, u_block,  offsets[1])

#  Initialize the state.
u0 = InitialCondition(num_equation)
sol = mfem.ParGridFunction(vfes, u_block.GetData())
sol.ProjectCoefficient(u0)

smyid = '{:0>6d}'.format(myid)
pmesh.Print("vortex-mesh."+smyid, 8)
for k in range(num_equation):
    uk = mfem.ParGridFunction(fes, u_block.GetBlock(k).GetData())
    sol_name = "vortex-" + str(k) + "-init."+smyid
    uk.Save(sol_name, 8)

# 9. Set up the nonlinear form corresponding to the DG discretization of the
#    flux divergence, and assemble the corresponding mass matrix.
Aflux = mfem.MixedBilinearForm(dfes, fes)
Aflux.AddDomainIntegrator(DomainIntegrator(dim))
Aflux.Assemble()

A = mfem.ParNonlinearForm(vfes)
rsolver = RiemannSolver()
ii = FaceIntegrator(rsolver, dim)
A.AddInteriorFaceIntegrator(ii)

# 10. Define the time-dependent evolution operator describing the ODE
#     right-hand side, and perform time-integration (looping over the time
#     iterations, ti, with a time-step dt).
euler = FE_Evolution(vfes, A, Aflux.SpMat())

if (visualization):
    MPI.COMM_WORLD.Barrier()
    sout = mfem.socketstream("localhost", 19916)
    sout.send_text("parallel " + str(num_procs) + " " + str(myid))
    sout.precision(8)
    sout.send_solution(pmesh, mom)
    sout.send_text("pause")
    sout.flush()
    if myid == 0:
        print("GLVis visualization paused.")
        print(" Press space (in the GLVis window) to resume it.")

# Determine the minimum element size.
my_hmin = 0
if (cfl > 0):
    my_hmin = min([pmesh.GetElementSize(i, 1) for i in range(pmesh.GetNE())])
hmin = MPI.COMM_WORLD.allreduce(my_hmin, op=MPI.MIN)

t = 0.0
euler.SetTime(t)
ode_solver.Init(euler)
if (cfl > 0):
    #  Find a safe dt, using a temporary vector. Calling Mult() computes the
    #  maximum char speed at all quadrature points on all faces.
    z = mfem.Vector(A.Width())
    A.Mult(sol, z)
    max_char_speed = MPI.COMM_WORLD.allreduce(
        ex18_common.max_char_speed, op=MPI.MAX)
    ex18_common.max_char_speed = max_char_speed
    dt = cfl * hmin / ex18_common.max_char_speed / (2*order+1)

# Integrate in time.
done = False
ti = 0
while not done:
    dt_real = min(dt, t_final - t)
    t, dt_real = ode_solver.Step(sol, t, dt_real)

    if (cfl > 0):
        max_char_speed = MPI.COMM_WORLD.allreduce(
            ex18_common.max_char_speed, op=MPI.MAX)
        ex18_common.max_char_speed = max_char_speed
        dt = cfl * hmin / ex18_common.max_char_speed / (2*order+1)

    ti = ti+1
    done = (t >= t_final - 1e-8*dt)
    if (done or ti % vis_steps == 0):
        if myid == 0:
            print("time step: " + str(ti) + ", time: " + "{:g}".format(t))
        if (visualization):
            sout.send_text("parallel " + str(num_procs) + " " + str(myid))
            sout.send_solution(pmesh, mom)
            sout.flush()

if myid == 0:
    print("done")

# 11. Save the final solution. This output can be viewed later using GLVis:
#     "glvis -np 4 -m vortex-mesh -g vortex-1-final".

for k in range(num_equation):
    uk = mfem.ParGridFunction(fes, u_block.GetBlock(k).GetData())
    sol_name = "vortex-" + str(k) + "-final."+smyid
    uk.Save(sol_name, 8)

# 12. Compute the L2 solution error summed for all components.
# if (t_final == 2.0):
if True:
    error = sol.ComputeLpError(2., u0)
    if myid == 0:
        print("Solution error: " + "{:g}".format(error))
