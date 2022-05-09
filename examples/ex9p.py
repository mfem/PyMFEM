'''
   MFEM example 9
      This is a version of Example 1 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname, exists
from mpi4py import MPI
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

num_proc = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
verbose = (myid == 0)

problem = 0
ser_ref_levels = 2
par_ref_levels = 0
order = 3
ode_solver_type = 4
t_final = 10
dt = 0.01
vis_steps = 5

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

# 3. Read the serial mesh from the given mesh file on all processors. We can
#    handle geometrically periodic meshes in this code.
meshfile = expanduser(join(path, 'data', 'periodic-hexagon.mesh'))
if not exists(meshfile):
    path = dirname(dirname(__file__))
    meshfile = expanduser(join(path, 'data', 'periodic-hexagon.mesh'))

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

# 5. Refine the mesh to increase the resolution. In this example we do
#    'ref_levels' of uniform refinement, where 'ref_levels' is a
#    command-line parameter. If the mesh is of NURBS type, we convert it to
#    a (piecewise-polynomial) high-order mesh.
for lev in range(ser_ref_levels):
    mesh.UniformRefinement()
    if mesh.NURBSext:
        mesh.SetCurvature(max(order, 1))
    bb_min, bb_max = mesh.GetBoundingBox(max(order, 1))


# 6. Define the parallel mesh by a partitioning of the serial mesh. Refine
#    this mesh further in parallel to increase the resolution. Once the
#    parallel mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
for k in range(par_ref_levels):
    pmesh.UniformRefinement()

# 7. Define the discontinuous DG finite element space of the given
#    polynomial order on the refined mesh.
fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
fes = mfem.ParFiniteElementSpace(pmesh, fec)

global_vSize = fes.GlobalTrueVSize()
if myid == 0:
    print("Number of unknowns: " + str(global_vSize))

#
#  Define coefficient using VecotrPyCoefficient and PyCoefficient
#  A user needs to define EvalValue method
#


class velocity_coeff(mfem.VectorPyCoefficient):
    def EvalValue(self, x):
        dim = len(x)

        center = (bb_min + bb_max)/2.0
        # map to the reference [-1,1] domain
        X = 2 * (x - center) / (bb_max - bb_min)
        if problem == 0:
            if dim == 1:
                v = [1.0, ]
            elif dim == 2:
                v = [sqrt(2./3.), sqrt(1./3)]
            elif dim == 3:
                v = [sqrt(3./6.), sqrt(2./6), sqrt(1./6.)]
        elif (problem == 1 or problem == 2):
            # Clockwise rotation in 2D around the origin
            w = pi/2
            if dim == 1:
                v = [1.0, ]
            elif dim == 2:
                v = [w*X[1],  - w*X[0]]
            elif dim == 3:
                v = [w*X[1],  - w*X[0],  0]
        elif (problem == 3):
            # Clockwise twisting rotation in 2D around the origin
            w = pi/2
            d = max((X[0]+1.)*(1.-X[0]), 0.) * max((X[1]+1.)*(1.-X[1]), 0.)
            d = d ** 2
            if dim == 1:
                v = [1.0, ]
            elif dim == 2:
                v = [d*w*X[1],  - d*w*X[0]]
            elif dim == 3:
                v = [d*w*X[1],  - d*w*X[0],  0]
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
                rx = 0.45
                ry = 0.25
                cx = 0.
                cy = -0.2
                w = 10.
                if dim == 3:
                    s = (1. + 0.25*cos(2 * pi * x[2]))
                    rx = rx * s
                    ry = ry * s
                return (erfc(w * (X[0]-cx-rx)) * erfc(-w*(X[0]-cx+rx)) *
                        erfc(w * (X[1]-cy-ry)) * erfc(-w*(X[1]-cy+ry)))/16

        elif problem == 2:
            rho = hypot(x[0], x[1])
            phi = arctan2(x[1], x[0])
            return (sin(pi * rho) ** 2) * sin(3*phi)
        elif problem == 3:
            return sin(pi * X[0]) * sin(pi * X[1])

        return 0.0

# Inflow boundary condition (zero for the problems considered in this example)


class inflow_coeff(mfem.PyCoefficient):
    def EvalValue(self, x):
        return 0

# 8. Set up and assemble the bilinear and linear forms corresponding to the
#    DG discretization. The DGTraceIntegrator involves integrals over mesh
#    interior faces.


velocity = velocity_coeff(dim)
inflow = inflow_coeff()
u0 = u0_coeff()

m = mfem.ParBilinearForm(fes)
m.AddDomainIntegrator(mfem.MassIntegrator())
k = mfem.ParBilinearForm(fes)
k.AddDomainIntegrator(mfem.ConvectionIntegrator(velocity, -1.0))
k.AddInteriorFaceIntegrator(
    mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5)))
k.AddBdrFaceIntegrator(
    mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5)))

b = mfem.ParLinearForm(fes)
b.AddBdrFaceIntegrator(
    mfem.BoundaryFlowIntegrator(inflow, velocity, -1.0, -0.5))

m.Assemble()
m.Finalize()
skip_zeros = 0
k.Assemble(skip_zeros)
k.Finalize(skip_zeros)
b.Assemble()

M = m.ParallelAssemble()
K = k.ParallelAssemble()
B = b.ParallelAssemble()

# 7. Define the initial conditions, save the corresponding grid function to
#    a file
u = mfem.ParGridFunction(fes)
u.ProjectCoefficient(u0)
U = u.GetTrueDofs()

smyid = '{:0>6d}'.format(myid)
mesh_name = "ex9-mesh."+smyid
sol_name = "ex9-init."+smyid
pmesh.Print(mesh_name, 8)
u.Save(sol_name, 8)


class FE_Evolution(mfem.PyTimeDependentOperator):
    def __init__(self, M, K, b):
        mfem.PyTimeDependentOperator.__init__(self, M.Height())

        self.M_prec = mfem.HypreSmoother()
        self.M_solver = mfem.CGSolver(M.GetComm())
        self.z = mfem.Vector(M.Height())

        self.K = K
        self.M = M
        self.b = b
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)
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
        self.K.Mult(x, self.z)
        self.z += b
        self.M_solver.Mult(self.z, y)


adv = FE_Evolution(M, K, B)

ode_solver.Init(adv)
t = 0.0
ti = 0
while True:
    if t > t_final - dt/2:
        break
    t, dt = ode_solver.Step(U, t, dt)
    ti = ti + 1
    if ti % vis_steps == 0:
        if myid == 0:
            print("time step: " + str(ti) + ", time: " + str(np.round(t, 3)))


u.Assign(U)
sol_name = "ex9-final."+smyid
u.Save(sol_name, 8)
