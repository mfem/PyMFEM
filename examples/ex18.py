'''
   MFEM example 18
      This is a version of Example 18 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

from .ex18_common import *

# Equation constant parameters.
num_equation = 4
specific_heat_ratio = 1.4
gas_constant = 1.0
        
parser = ArgParser(description='Ex18')
parser.add_argument('-m', '--mesh',
                    default = 'periodic-square.mesh',
                    action = 'store', type = str,
                    help='Mesh file to use.')
parser.add_argument('-p', '--problem',
                    action = 'store', default = 1, type=int,                    
                    help = 'Problem setup to use. See options in velocity_function().')
parser.add_argument('-r', '--refine',
                    action = 'store', default = 1, type=int,
                    help = "Number of times to refine the mesh uniformly.")
parser.add_argument('-o', '--order',
                    action = 'store', default = 1, type=int,
                    help = "Finite element order (polynomial degree)")
parser.add_argument('-s', '--ode_solver',
                    action = 'store', default = 4, type=int,                    
                    help = "ODE solver: 1 - Forward Euler,\n\t" + 
                           "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6.")
parser.add_argument('-t', '--t-final',
                    action = 'store', default = 2.0, type=float,
                    help = "Final time; start time is 0.")   
parser.add_argument("-dt", "--time-step",
                    action = 'store', default = -0.01, type=float,
                    help = "Time step.");
parser.add_argument('-c', '--cfl-number',
                    action = 'store', default = 0.3, type=float,
                    help = "CFL number for timestep calculation.")
parser.add_argument('-vis', '--visualization',
                    action = 'store_true',
                    help='Enable GLVis visualization')
parser.add_argument('-vs', '--visualization-steps',
                    action = 'store', default = 50, type=float,
                    help = "Visualize every n-th timestep.")

args = parser.parse_args()
mesh = args.mesh
ref_levels = args.refine
order = args.order
ode_solver_type = args.ode_solver
t_final = args.t_final
dt = args.time_step
cfl = args.cfl_number
visualization = args.visualization
vis_stesp = args.visualization_steps
problem = args.problem

# 2. Read the mesh from the given mesh file. This example requires a 2D
#    periodic mesh, such as ../data/periodic-square.mesh.
meshfile = expanduser(join(path, 'data', mesh))
mesh = mfem.Mesh(meshfile, 1,1)
dim = mesh.Dimension()

# 3. Define the ODE solver used for time integration. Several explicit
#    Runge-Kutta methods are available.
ode_solver = None
if ode_solver_type == 1:   ode_solver = mfem.ForwardEulerSolver()
elif ode_solver_type == 2: ode_solver = mfem.RK2Solver(1.0)
elif ode_solver_type == 3: ode_solver = mfem.RK3SSolver()
elif ode_solver_type == 4: ode_solver = mfem.RK4Solver()
elif ode_solver_type == 6: ode_solver = mfem.RK6Solver()
else:
    print( "Unknown ODE solver type: " + str(ode_solver_type))
    exit

# 4. Refine the mesh to increase the resolution. In this example we do
#    'ref_levels' of uniform refinement, where 'ref_levels' is a
#    command-line parameter.
for lev in range(ref_levels):
    mesh.UniformRefinement();

# 5. Define the discontinuous DG finite element space of the given
#    polynomial order on the refined mesh.

fec = mfem.DG_FECollection(order, dim)
# Finite element space for a scalar (thermodynamic quantity)
fes = mfem.FiniteElementSpace(mesh, fec)
# Finite element space for a mesh-dim vector quantity (momentum)
dfes = mfem.FiniteElementSpace(mesh, fec, dim, mfem.Ordering.byNODES)
# Finite element space for all variables together (total thermodynamic state)
vfes = mfem.FiniteElementSpace(mesh, fec, num_equation, mfem.Ordering.byNODES)

assert fes.GetOrdering() == mfem.Ordering.byNODES, "Ordering must be byNODES")
print("Number of unknowns: " + str(fes.GetVSize()))

# 6. Define the initial conditions, save the corresponding mesh and grid
#    functions to a file. This can be opened with GLVis with the -gc option.
#    The solution u has components {density, x-momentum, y-momentum, energy}.
#    These are stored contiguously in the BlockVector u_block.
offsets = [k*vfes.GetNDofs() for k in range(num_equation+1)]
offsets(num_equation + 1)
u_block = mfem.BlockVector(offsets)
mfem.Vector(u_block, 
mom = mfem.GridFunction mom(dfes, u_block,  offsets[1])
            
#
#  Define coefficient using VecotrPyCoefficient and PyCoefficient
#  A user needs to define EvalValue method        
#

u0 = mfem.VectorFunctionCoefficient(num_equation, InitialCondition)
sol = mfem.GridFunction(vfes, u_block.GetData())
sol.ProjectCoeffcient(u0)
            

#  7. Set up the nonlinear form corresponding to the DG discretization of the
#     flux divergence, and assemble the corresponding mass matrix.
Aflux = mfem.MixedBilinearForm(dfes, fes)
Aflux.AddDomainIntegrator(DomainIntegrator(dim))
Aflux.Assemble();

A = mfem.NonlinearForm(vfes)
rsolver = RiemannSolver()
A.AddInteriorFaceIntegrator(FaceIntegrator(rsolver, dim))

#  8. Define the time-dependent evolution operator describing the ODE
#     right-hand side, and perform time-integration (looping over the time
#     iterations, ti, with a time-step dt).
euler = FE_Evolution(vfes, A, Aflux.SpMat())
            
# Determine the minimum element size.
double hmin;
if (cfl > 0):
   hmin = min([mesh.GetElementSize(i, 1) for i in range(mesh.GetNE())])

t = 0.0            
euler.SetTime(t);
ode_solver.Init(euler);
if (cfl > 0):
    #  Find a safe dt, using a temporary vector. Calling Mult() computes the
    #  maximum char speed at all quadrature points on all faces.
    z = mfem.Vector(A.Width())
    max_char_speed = 0.
    z =  A.Mult(sol)
    dt = cfl * hmin / max_char_speed / (2*order+1)

# Integrate in time.
done = false;
ti = 0            
while not done:
    dt_real = min(dt, t_final - t);
    ode_solver.Step(sol, t, dt_real);
    if (cfl > 0):
         dt = cfl * hmin / max_char_speed / (2*order+1);
    ti++
    done = (t >= t_final - 1e-8*dt)
    if (done || ti % vis_steps == 0):
        print("time step: " + str(ti) + ", time: " + str(t))
        if (visualization):
            sout << "solution\n" << mesh << mom << flush;
                        
class velocity_coeff(mfem.VectorPyCoefficient):
   def EvalValue(self, x):        
       dim = len(x)
        
       center = (bb_min + bb_max)/2.0
       # map to the reference [-1,1] domain                
       X = 2 * (x - center) / (bb_max - bb_min)
       if problem == 0:
           if dim == 1: v = [1.0,]
           elif dim == 2: v = [sqrt(2./3.), sqrt(1./3)]
           elif dim == 3: v = [sqrt(3./6.), sqrt(2./6), sqrt(1./6.)]
       elif (problem == 1 or problem == 2):
           # Clockwise rotation in 2D around the origin                
           w = pi/2
           if dim == 1: v = [1.0,]
           elif dim == 2: v = [w*X[1],  - w*X[0]]
           elif dim == 3: v = [w*X[1],  - w*X[0],  0]
       elif (problem == 3):
           # Clockwise twisting rotation in 2D around the origin
           w = pi/2
           d = max((X[0]+1.)*(1.-X[0]),0.) * max((X[1]+1.)*(1.-X[1]),0.)
           d = d ** 2
           if dim == 1: v = [1.0,]
           elif dim == 2: v = [d*w*X[1],  - d*w*X[0]]
           elif dim == 3: v = [d*w*X[1],  - d*w*X[0],  0]
       return v
        
class u0_coeff(mfem.PyCoefficient):
   def EvalValue(self, x):        
       dim = len(x)

       center = (bb_min + bb_max)/2.0
       # map to the reference [-1,1] domain        
       X = 2 * (x - center) / (bb_max - bb_min)
       if (problem == 0 or problem == 1):
           if dim == 1:
               return exp(-40. * (X[0]-0.5)**2)
           elif (dim == 2 or dim == 3):
               rx = 0.45; ry = 0.25; cx = 0.; cy = -0.2; w = 10.
               if dim == 3:
                   s = (1. + 0.25*cos(2 * pi * x[2]))
                   rx = rx * s
                   ry = ry * s
               return ( erfc( w * (X[0]-cx-rx)) * erfc(-w*(X[0]-cx+rx)) *
                        erfc( w * (X[1]-cy-ry)) * erfc(-w*(X[1]-cy+ry)) )/16
        
       elif problem == 2:
           rho = hypot(x[0], x[1])
           phi = arctan2(x[1], x[0])
           return (sin(pi * rho) **2) * sin(3*phi)
       elif problem == 3:
           return sin(pi * X[0]) * sin(pi * X[1])
       
       return 0.0

# Inflow boundary condition (zero for the problems considered in this example)
class inflow_coeff(mfem.PyCoefficient):
   def EvalValue(self, x):
       return 0

# 6. Set up and assemble the bilinear and linear forms corresponding to the
#    DG discretization. The DGTraceIntegrator involves integrals over mesh
#    interior faces.

velocity = velocity_coeff(dim)
inflow = inflow_coeff()
u0 = u0_coeff()

m = mfem.BilinearForm(fes)
m.AddDomainIntegrator(mfem.MassIntegrator())
k = mfem.BilinearForm(fes)
k.AddDomainIntegrator(mfem.ConvectionIntegrator(velocity, -1.0))
k.AddInteriorFaceIntegrator(
      mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5)))
k.AddBdrFaceIntegrator(
      mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5)))

b = mfem.LinearForm(fes)
b.AddBdrFaceIntegrator(
      mfem.BoundaryFlowIntegrator(inflow, velocity, -1.0, -0.5))

m.Assemble()
m.Finalize()
skip_zeros = 0
k.Assemble(skip_zeros)
k.Finalize(skip_zeros)
b.Assemble()

# 7. Define the initial conditions, save the corresponding grid function to
#    a file 
u = mfem.GridFunction(fes)
u.ProjectCoefficient(u0)

mesh.PrintToFile('ex9.mesh', 8)
u.SaveToFile('ex9-init.gf', 8)

class FE_Evolution(mfem.PyTimeDependentOperator):
    def __init__(self, M, K, b):
        mfem.PyTimeDependentOperator.__init__(self, M.Size())        

        
        self.K = K        
        self.M = M
        self.b = b
        self.z = mfem.Vector(M.Size())
        self.zp = np.zeros(M.Size())
        self.M_prec = mfem.DSmoother()        
        self.M_solver = mfem.CGSolver()
        self.M_solver.SetPreconditioner(self.M_prec)        
        self.M_solver.SetOperator(M)
        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(1e-9)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(100)
        self.M_solver.SetPrintLevel(0)
        

#    def EvalMult(self, x):
#        if you want to impolement Mult in using python objects, 
#        such as numpy.. this needs to be implemented and don't
#        overwrite Mult

    def Mult(self, x, y):
        self.K.Mult(x, self.z);
        self.z += b;
        self.M_solver.Mult(self.z, y)

adv = FE_Evolution(m.SpMat(), k.SpMat(), b)

ode_solver.Init(adv)
t = 0.0; ti = 0;
while True:
   if t > t_final - dt/2: break
   t, dt = ode_solver.Step(u, t, dt);
   ti = ti + 1
   if ti % vis_steps == 0:
       print("time step: " + str(ti) + ", time: " + str(t))

u.SaveToFile('ex9-final.gf', 8)
