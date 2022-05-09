'''
   MFEM example 10

      This examples solves a time dependent nonlinear elasticity
      problem of the form dv/dt = H(x) + S v, dx/dt = v, where H is a
      hyperelastic model and S is a viscosity operator of Laplacian
      type.

      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
import sys
from mfem.common.arg_parser import ArgParser
import mfem.ser as mfem
from mfem.ser import intArray, add_vector, Add
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

parser = ArgParser(description='Ex10')
parser.add_argument('-m', '--mesh',
                    default='beam-quad.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-r', '--refine-serial',
                    action='store', default=2, type=int,
                    help="Number of times to refine the mesh uniformly before parallel")
parser.add_argument('-o', '--order',
                    action='store', default=2, type=int,
                    help="Finite element order (polynomial degree)")
help_ode = "\n".join(["ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3",
                      "\t11 - Forward Euler, 12 - RK2",
                      "\t13 - RK3 SSP, 14 - RK4."])
parser.add_argument('-s', '--ode-solver',
                    action='store', default=3, type=int,
                    help=help_ode)
parser.add_argument('-tf', '--t-final',
                    action='store', default=300.0, type=float,
                    help="Final time; start time is 0.")
parser.add_argument('-dt', '--time-step',
                    action='store', default=3.0, type=float,
                    help="Time step")
parser.add_argument("-v", "--viscosity",
                    action='store', default=1e-2, type=float,
                    help="Viscosity coefficient.")
parser.add_argument("-mu", "--shear-modulus",
                    action='store', default=0.25, type=float,
                    help="Shear modulus in the Neo-Hookean hyperelastic model.")
parser.add_argument("-K", "--bulk-modulus",
                    action='store', default=5.0, type=float,
                    help="Bulk modulus in the Neo-Hookean hyperelastic model.")
parser.add_argument('-vis', '--visualization',
                    action='store_true', default=True,
                    help='Enable GLVis visualization')
parser.add_argument("-vs", "--visualization-steps",
                    action='store', default=1,  type=int,
                    help="Visualize every n-th timestep.")
args = parser.parse_args()

ref_levels = args.refine_serial
order = args.order
ode_solver_type = args.ode_solver
t_final = args.t_final
dt = args.time_step
visc = args.viscosity
mu = args.shear_modulus
K = args.bulk_modulus
visualization = args.visualization
vis_steps = args.visualization_steps
meshfile = expanduser(
    join(dirname(__file__), '..', 'data', args.mesh))

parser.print_options(args)

mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
#        self.solver.SetOperator(M)
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

for lev in range(ref_levels):
    mesh.UniformRefinement()

# 5. Define the vector finite element spaces representing the mesh
#    deformation x, the velocity v, and the initial configuration, x_ref.
#    Define also the elastic energy density, w, which is in a discontinuous
#    higher-order space. Since x and v are integrated in time as a system,
#    we group them together in block vector vx, with offsets given by the
#    fe_offset array.
fec = mfem.H1_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fec, dim)

fe_size = fespace.GetVSize()
print("Number of velocity/deformation unknowns: " + str(fe_size))
fe_offset = intArray([0, fe_size, 2*fe_size])

vx = mfem.BlockVector(fe_offset)
x = mfem.GridFunction()
v = mfem.GridFunction()
v.MakeRef(fespace, vx.GetBlock(0), 0)
x.MakeRef(fespace, vx.GetBlock(1), 0)

x_ref = mfem.GridFunction(fespace)
mesh.GetNodes(x_ref)

w_fec = mfem.L2_FECollection(order + 1, dim)
w_fespace = mfem.FiniteElementSpace(mesh, w_fec)
w = mfem.GridFunction(w_fespace)

# 6. Set the initial conditions for v and x, and the boundary conditions on
#    a beam-like mesh (see description above).


class InitialVelocity(mfem.VectorPyCoefficient):
    def EvalValue(self, x):
        dim = len(x)
        s = 0.1/64.

        v = np.zeros(len(x))
        v[-1] = s*x[0]**2*(8.0-x[0])
        v[0] = -s*x[0]**2
        return v


class InitialDeformation(mfem.VectorPyCoefficient):
    def EvalValue(self, x):
        return x.copy()


velo = InitialVelocity(dim)
v.ProjectCoefficient(velo)
deform = InitialDeformation(dim)
x.ProjectCoefficient(deform)

ess_bdr = intArray(fespace.GetMesh().bdr_attributes.Max())
ess_bdr.Assign(0)
ess_bdr[0] = 1

# 7. Define HyperelasticOperator and initialize it
#    the initial energies.


class ElasticEnergyCoefficient(mfem.PyCoefficient):
    def __init__(self, model, x):
        self.x = x
        self.model = model
        self.J = mfem.DenseMatrix()
        mfem.PyCoefficient.__init__(self)

    def Eval(self, T, ip):
        self.model.SetTransformation(T)
        self.x.GetVectorGradient(T, self.J)
        # T.Jacobian().Print()
        # print self.x.GetDataArray()
        # self.J.Print()
        return self.model.EvalW(self.J)/(self.J.Det())


class ReducedSystemOperator(mfem.PyOperator):
    def __init__(self, M, S, H):
        mfem.PyOperator.__init__(self, M.Height())
        self.M = M
        self.S = S
        self.H = H
        self.Jacobian = None
        h = M.Height()
        self.w = mfem.Vector(h)
        self.z = mfem.Vector(h)
        self.dt = 0.0
        self.v = None
        self.x = None

    def SetParameters(self, dt, v, x):
        self.dt = dt
        self.v = v
        self.x = x

    def Mult(self, k, y):
        add_vector(self.v, self.dt, k, self.w)
        add_vector(self.x, self.dt, self.w, self.z)
        self.H.Mult(self.z, y)
        self.M.AddMult(k, y)
        self.S.AddMult(self.w, y)

    def GetGradient(self, k):
        Jacobian = Add(1.0, self.M.SpMat(), self.dt, self.S.SpMat())
        self.Jacobian = Jacobian
        add_vector(self.v, self.dt, k, self.w)
        add_vector(self.x, self.dt, self.w, self.z)
        grad_H = self.H.GetGradientMatrix(self.z)

        Jacobian.Add(self.dt**2, grad_H)
        return Jacobian


class HyperelasticOperator(mfem.PyTimeDependentOperator):
    def __init__(self, fespace, ess_bdr, visc, mu, K):
        mfem.PyTimeDependentOperator.__init__(self, 2*fespace.GetVSize(), 0.0)

        rel_tol = 1e-8
        skip_zero_entries = 0
        ref_density = 1.0
        self.z = mfem.Vector(self.Height()//2)
        self.fespace = fespace
        self.viscosity = visc

        M = mfem.BilinearForm(fespace)
        S = mfem.BilinearForm(fespace)
        H = mfem.NonlinearForm(fespace)
        self.M = M
        self.H = H
        self.S = S

        rho = mfem.ConstantCoefficient(ref_density)
        M.AddDomainIntegrator(mfem.VectorMassIntegrator(rho))
        M.Assemble(skip_zero_entries)
        M.EliminateEssentialBC(ess_bdr)
        M.Finalize(skip_zero_entries)

        M_solver = mfem.CGSolver()
        M_prec = mfem.DSmoother()
        M_solver.iterative_mode = False
        M_solver.SetRelTol(rel_tol)
        M_solver.SetAbsTol(0.0)
        M_solver.SetMaxIter(30)
        M_solver.SetPrintLevel(0)
        M_solver.SetPreconditioner(M_prec)
        M_solver.SetOperator(M.SpMat())

        self.M_solver = M_solver
        self.M_prec = M_prec

        model = mfem.NeoHookeanModel(mu, K)
        H.AddDomainIntegrator(mfem.HyperelasticNLFIntegrator(model))
        H.SetEssentialBC(ess_bdr)
        self.model = model

        visc_coeff = mfem.ConstantCoefficient(visc)
        S.AddDomainIntegrator(mfem.VectorDiffusionIntegrator(visc_coeff))
        S.Assemble(skip_zero_entries)
        S.EliminateEssentialBC(ess_bdr)
        S.Finalize(skip_zero_entries)

        self.reduced_oper = ReducedSystemOperator(M, S, H)

        J_prec = mfem.DSmoother(1)
        J_minres = mfem.MINRESSolver()
        J_minres.SetRelTol(rel_tol)
        J_minres.SetAbsTol(0.0)
        J_minres.SetMaxIter(300)
        J_minres.SetPrintLevel(-1)
        J_minres.SetPreconditioner(J_prec)

        self.J_solver = J_minres
        self.J_prec = J_prec

        newton_solver = mfem.NewtonSolver()
        newton_solver.iterative_mode = False
        newton_solver.SetSolver(self.J_solver)
        newton_solver.SetOperator(self.reduced_oper)
        newton_solver.SetPrintLevel(1)  # print Newton iterations
        newton_solver.SetRelTol(rel_tol)
        newton_solver.SetAbsTol(0.0)
        newton_solver.SetMaxIter(10)
        self.newton_solver = newton_solver

    def Mult(self, vx, vx_dt):
        sc = self.Height()//2
        v = mfem.Vector(vx, 0,  sc)
        x = mfem.Vector(vx, sc,  sc)
        dv_dt = mfem.Vector(dvx_dt, 0, sc)
        dx_dt = mfem.Vector(dvx_dt, sc,  sc)
        self.H.Mult(x, z)
        if (self.viscosity != 0.0):
            S.AddMult(v, z)
        z.Neg()
        M_solver.Mult(z, dv_dt)
        dx_dt = v
#        Print(vx.Size())

    def ImplicitSolve(self, dt, vx, dvx_dt):
        sc = self.Height()//2
        v = mfem.Vector(vx, 0,  sc)
        x = mfem.Vector(vx, sc,  sc)
        dv_dt = mfem.Vector(dvx_dt, 0, sc)
        dx_dt = mfem.Vector(dvx_dt, sc,  sc)

        # By eliminating kx from the coupled system:
        # kv = -M^{-1}*[H(x + dt*kx) + S*(v + dt*kv)]
        # kx = v + dt*kv
        # we reduce it to a nonlinear equation for kv, represented by the
        # backward_euler_oper. This equation is solved with the newton_solver
        # object (using J_solver and J_prec internally).
        self.reduced_oper.SetParameters(dt, v, x)
        zero = mfem.Vector()  # empty vector is interpreted as
        # zero r.h.s. by NewtonSolver
        self.newton_solver.Mult(zero, dv_dt)
        add_vector(v, dt, dv_dt, dx_dt)

    def ElasticEnergy(self, x):
        return self.H.GetEnergy(x)

    def KineticEnergy(self, v):
        return 0.5*self.M.InnerProduct(v, v)

    def GetElasticEnergyDensity(self, x, w):
        w_coeff = ElasticEnergyCoefficient(self.model, x)
        w.ProjectCoefficient(w_coeff)


oper = HyperelasticOperator(fespace, ess_bdr, visc, mu, K)
ee0 = oper.ElasticEnergy(x)
ke0 = oper.KineticEnergy(v)

print("initial elastic energy (EE) = " + str(ee0))
print("initial kinetic energy (KE) = " + str(ke0))
print("initial   total energy (TE) = " + str(ee0 + ke0))

# 8. Perform time-integration (looping over the time iterations, ti, with a
#    time-step dt).
ode_solver.Init(oper)
t = 0.
ti = 1
last_step = False
while not last_step:
    if (t + dt >= t_final - dt/2):
        last_step = True

    t, dt = ode_solver.Step(vx, t, dt)

    if (last_step or (ti % vis_steps) == 0):
        ee = oper.ElasticEnergy(x)
        ke = oper.KineticEnergy(v)

        text = ("step " + str(ti) + ", t = " + str(t) + ", EE = " +
                str(ee) + ", KE = " + str(ke) +
                ", dTE = " + str((ee+ke)-(ee0+ke0)))

        print(text)
    ti = ti + 1
#
# if i translate c++ line-by-line, ti seems the second swap does not work...
#
nodes = x
owns_nodes = 0
nodes, owns_nodes = mesh.SwapNodes(nodes, owns_nodes)
mesh.Print('deformed.mesh', 8)
mesh.SwapNodes(nodes, owns_nodes)
v.Save('velocity.sol', 8)
oper.GetElasticEnergyDensity(x, w)
w.Save('elastic_energy.sol',  8)
