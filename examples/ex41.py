'''
   MFEM example 41 (converted from ex41.cpp)

   See c++ version in the MFEM library for more detail

   Sample runs:
      python ex41.py
      python ex41.py -cg
      python ex41.py -m ../data/periodic-hexagon.mesh -p 0 -r 2 -dt 0.005 -tf 10
      python ex41.py -m ../data/periodic-square.mesh -p 1 -r 2 -dt 0.005 -tf 9
      python ex41.py -m ../data/periodic-hexagon.mesh -p 1 -r 2 -dt 0.005 -tf 9
      python ex41.py -m ../data/amr-quad.mesh -p 1 -r 2 -dt 0.002 -tf 9
      python ex41.py -m ../data/star-q3.mesh -p 1 -r 2 -dt 0.001 -tf 9
      python ex41.py -m ../data/star-mixed.mesh -p 1 -r 2 -dt 0.005 -tf 9
      python ex41.py -m ../data/disc-nurbs.mesh -p 1 -r 3 -dt 0.005 -tf 9
      python ex41.py -m ../data/disc-nurbs.mesh -p 2 -r 3 -dt 0.005 -tf 9
      python ex41.py -m ../data/periodic-square.mesh -p 3 -r 4 -dt 0.0025 -tf 9 -vs 20
      python ex41.py -m ../data/periodic-cube.mesh -p 0 -r 2 -o 2 -dt 0.01 -tf 8

   Description:  This example code solves the time-dependent advection-diffusion
                 equation du/dt + v.grad(u) - a div(grad(u)) = 0, where v is a
                 given fluid velocity, a is the diffusion coefficient, and
                 u0(x)=u(0,x) is a given initial condition.

                 The example demonstrates the use of Discontinuous Galerkin (DG)
                 bilinear forms in MFEM (face integrators), and the use of IMEX
                 ODE time integrators.

                 The option to use continuous finite elements is available too.


   Python note: binary (Sidre) VisIt output is not translated.
'''
import os
from os.path import expanduser, join
from math import erfc
import numpy as np

import mfem.ser as mfem


class Implicit_Solver(mfem.Solver):
    """Solver for the implicit part of the ODE (the diffusion term).

    Solves systems of the form: (M + dt*S) k = rhs.
    """
    def __init__(self, M, S, fes):
        mfem.Solver.__init__(self, M.Height())
        self.M = M
        self.S = S
        self.A = None
        self.dt = None
        self.prec = mfem.BlockILU(
            fes.GetTypicalFE().GetDof(),
            mfem.BlockILU.Reordering_MINIMUM_DISCARDED_FILL)
        self.linear_solver = mfem.CGSolver()
        self.linear_solver.iterative_mode = False
        self.linear_solver.SetRelTol(1e-9)
        self.linear_solver.SetAbsTol(0.0)
        self.linear_solver.SetMaxIter(100)
        self.linear_solver.SetPrintLevel(0)
        self.linear_solver.SetPreconditioner(self.prec)

    def SetTimeStep(self, dt):
        if self.dt is None or abs(self.dt - dt) > 10*np.finfo(float).eps:
            self.dt = dt
            # Form operator A = M + dt*S
            self.A = mfem.Add(1.0, self.M, dt, self.S)
            # This will also call SetOperator on the preconditioner.
            self.linear_solver.SetOperator(self.A)

    def SetOperator(self, op):
        self.linear_solver.SetOperator(op)

    def Mult(self, x, y):
        self.linear_solver.Mult(x, y)


class IMEX_Evolution(mfem.PyTimeDependentOperator):
    """A time-dependent operator for the right-hand side of the ODE.

    The weak form of the advection-diffusion equation is M du/dt = K u - S u + b,
    where M is the mass matrix, K and S are the advection and diffusion
    matrices, and b describes the flow on the boundary. In the case of IMEX
    evolution, the diffusion term is treated implicitly, and the advection
    term is treated explicitly.
    """
    def __init__(self, M, K, S, b):
        mfem.PyTimeDependentOperator.__init__(self, M.FESpace().GetTrueVSize())
        self.M, self.K, self.S, self.b = M, K, S, b
        self.z = mfem.Vector(self.Height())
        self.M_prec = mfem.DSmoother(M.SpMat())
        self.M_solver = mfem.CGSolver()
        self.M_solver.SetOperator(M.SpMat())
        self.implicit_solver = Implicit_Solver(M.SpMat(), S.SpMat(), M.FESpace())
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(1e-9)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(100)
        self.M_solver.SetPrintLevel(0)

    def Mult1(self, x, y):
        # Perform the explicit step
        # y = M^{-1} (K x + b)
        self.K.Mult(x, self.z)
        self.z += self.b
        self.M_solver.Mult(self.z, y)

    def ImplicitSolve2(self, dt, x, k):
        # Perform the implicit step
        # solve for k, k = -(M+dt S)^{-1} S x
        self.S.Mult(x, self.z)
        self.z.Neg()
        self.implicit_solver.SetTimeStep(dt)
        self.implicit_solver.Mult(self.z, k)

    def Mult(self, x, y):
        if self.GetEvalMode() == mfem.TimeDependentOperator.ADDITIVE_TERM_1:
            self.Mult1(x, y)
        else:
            raise RuntimeError('TimeDependentOperator::Mult() is not overridden!')

    def ImplicitSolve(self, dt, x, k):
        if self.GetEvalMode() == mfem.TimeDependentOperator.ADDITIVE_TERM_2:
            self.ImplicitSolve2(dt, x, k)
        else:
            raise RuntimeError('TimeDependentOperator::ImplicitSolve() is not overridden!')


def run(meshfile='', problem=0, ref_levels=2, order=3,
        ode_solver_type=64, t_final=10.0, dt=0.01, diffusion_term=0.01,
        cg=False, vis_steps=50, visualization=True, visit=False, paraview=False):

    # 2. Read the mesh from the given mesh file. We can handle geometrically
    #    periodic meshes in this code.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    # 3. Define the IMEX (Split) ODE solver used for time integration. The IMEX
    #    solvers currently available are: 61 - Forward Backward Euler,
    #    62 - IMEXRK2(2,2,2), 63 - IMEXRK2(2,3,2), and 64 - IMEX_DIRK_RK3.
    if ode_solver_type == 61:
        ode_solver = mfem.IMEXExpImplEuler()
    elif ode_solver_type == 62:
        ode_solver = mfem.IMEXRK2()
    elif ode_solver_type == 63:
        ode_solver = mfem.IMEXRK2_3StageExplicit()
    elif ode_solver_type == 64:
        ode_solver = mfem.IMEX_DIRK_RK3()
    else:
        raise ValueError('Unknown IMEX ODE solver type: ' + str(ode_solver_type))

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement, where 'ref_levels' is a
    #    command-line parameter.
    for lev in range(ref_levels):
        mesh.UniformRefinement()
    if mesh.NURBSext:
        mesh.SetCurvature(max(order, 1))
    bb_min, bb_max = mesh.GetBoundingBox(max(order, 1))

    # 5. Define the discontinuous DG finite element space of the given
    #    polynomial order on the refined mesh.
    if cg:
        fec = mfem.H1_FECollection(order, dim)
    else:
        fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
    fes = mfem.FiniteElementSpace(mesh, fec)
    print('Number of unknowns: ' + str(fes.GetVSize()))

    # Velocity coefficient
    @mfem.jit.vector(shape=(dim,))
    def velocity(x):
        # Map to the reference [-1,1] domain.
        center = (bb_min + bb_max)*0.5
        X = 2*(x - center)/(bb_max - bb_min)
        v = np.zeros(dim)
        if problem == 0:
            # Translations in 1D, 2D, and 3D
            if dim == 1:
                v[0] = 1.0
            elif dim == 2:
                v[:] = [np.sqrt(2./3.), np.sqrt(1./3.)]
            else:
                v[:] = [np.sqrt(3./6.), np.sqrt(2./6.), np.sqrt(1./6.)]
        else:
            # Clockwise rotation in 2D around the origin
            w = np.pi/2
            if dim == 1:
                v[0] = 1.0
            else:
                d = 1.0
                if problem == 3:
                    # Clockwise twisting rotation in 2D around the origin
                    d = (max((X[0]+1)*(1-X[0]), 0.) *
                         max((X[1]+1)*(1-X[1]), 0.))**2
                v[0], v[1] = d*w*X[1], -d*w*X[0]
        return v

    # Initial condition
    @mfem.jit.scalar
    def u0(x):
        # Map to the reference [-1,1] domain.
        center = (bb_min + bb_max)*0.5
        X = 2*(x - center)/(bb_max - bb_min)
        if problem in (0, 1):
            if dim == 1:
                return np.exp(-40*(X[0]-0.5)**2)
            rx, ry, cx, cy, w = 0.45, 0.25, 0., -0.2, 10.
            if dim == 3:
                scale = 1 + 0.25*np.cos(2*np.pi*X[2])
                rx *= scale
                ry *= scale
            return (erfc(w*(X[0]-cx-rx))*erfc(-w*(X[0]-cx+rx)) *
                    erfc(w*(X[1]-cy-ry))*erfc(-w*(X[1]-cy+ry)))/16
        if problem == 2:
            rho = np.hypot(X[0], X[1])
            phi = np.arctan2(X[1], X[0])
            return np.sin(np.pi*rho)**2*np.sin(3*phi)
        return np.sin(np.pi*X[0])*np.sin(np.pi*X[1])

    # 6. Set up and assemble the bilinear and linear forms corresponding to the
    #    DG discretization. The DGTraceIntegrator involves integrals over mesh
    #    interior faces.
    diff_coeff = mfem.ConstantCoefficient(diffusion_term)
    m, k, s = mfem.BilinearForm(fes), mfem.BilinearForm(fes), mfem.BilinearForm(fes)
    b = mfem.Vector(fes.GetTrueVSize())
    b.Assign(0.0)  # The inflow on the boundaries is set to zero.
    m.AddDomainIntegrator(mfem.MassIntegrator())
    alpha, sigma, kappa = -1.0, -1.0, float((order+1)**2)
    k.AddDomainIntegrator(mfem.ConvectionIntegrator(velocity, alpha))
    s.AddDomainIntegrator(mfem.DiffusionIntegrator(diff_coeff))
    if not cg:
        k.AddInteriorFaceIntegrator(mfem.NonconservativeDGTraceIntegrator(velocity, alpha))
        k.AddBdrFaceIntegrator(mfem.NonconservativeDGTraceIntegrator(velocity, alpha))
        s.AddInteriorFaceIntegrator(mfem.DGDiffusionIntegrator(diff_coeff, sigma, kappa))
        s.AddBdrFaceIntegrator(mfem.DGDiffusionIntegrator(diff_coeff, sigma, kappa))
    skip_zeros = 0
    for form in (m, k, s):
        form.Assemble(skip_zeros)
        form.Finalize(skip_zeros)

    # 7. Define the initial conditions.
    u = mfem.GridFunction(fes)
    u.ProjectCoefficient(u0)

    # Create data collection for solution output.
    # Python note: only ASCII VisIt output is provided here.
    precision = 8
    if visit:
        dc = mfem.VisItDataCollection('Example41', mesh)
        dc.SetPrecision(precision)
        dc.RegisterField('solution', u)
        dc.SetCycle(0)
        dc.SetTime(0.0)
        dc.Save()

    # 8. Set up paraview visualization, if desired.
    if paraview:
        pv = mfem.ParaViewDataCollection('Example41', mesh)
        pv.SetPrefixPath('ParaView')
        pv.RegisterField('solution', u)
        pv.SetLevelsOfDetail(order)
        pv.SetDataFormat(mfem.VTKFormat_BINARY)
        pv.SetHighOrderOutput(True)
        pv.SetCycle(0)
        pv.SetTime(0.0)
        pv.Save()

    if visualization:
        sout = mfem.socketstream('localhost', 19916)
        if not sout.good():
            print('Unable to connect to GLVis server at localhost:19916')
            visualization = False
            print('GLVis visualization disabled.')
        else:
            sout.precision(precision)
            sout << 'solution\n' << mesh << u << 'pause\n'
            sout.flush()
            print('GLVis visualization paused. Press space (in the GLVis window) to resume it.')

    # 9. Define the time-dependent evolution operator describing the ODE
    #    right-hand side, and perform time-integration (looping over the time
    #    iterations, ti, with a time-step dt).
    adv = IMEX_Evolution(m, k, s, b)
    t = 0.0
    adv.SetTime(t)
    ode_solver.Init(adv)
    ti = 0
    while t < t_final - 1e-8*dt:
        dt_real = min(dt, t_final - t)
        t, dt_real = ode_solver.Step(u, t, dt_real)
        ti += 1
        done = t >= t_final - 1e-8*dt
        if done or ti % vis_steps == 0:
            print(f'time step: {ti}, time: {t:.6g}')
            if paraview:
                pv.SetCycle(ti)
                pv.SetTime(t)
                pv.Save()
            if visualization:
                sout << 'solution\n' << mesh << u
                sout.flush()
            if visit:
                dc.SetCycle(ti)
                dc.SetTime(t)
                dc.Save()


if __name__ == '__main__':
    from mfem.common.arg_parser import ArgParser

    # 1. Parse command-line options.
    parser = ArgParser(description='Ex41 (IMEX advection-diffusion)')
    parser.add_argument('-m', '--mesh', default='periodic-square.mesh', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-p', '--problem', default=0, type=int, choices=range(4),
                        help='Problem setup to use. See velocity_coeff.')
    parser.add_argument('-r', '--refine', default=2, type=int,
                        help='Number of times to refine the mesh uniformly.')
    parser.add_argument('-o', '--order', default=3, type=int,
                        help='Order of the finite elements.')
    parser.add_argument('-s', '--ode-solver', default=64, type=int,
                        help='61 - Forward Backward Euler, 62 - IMEXRK2(2,2,2), '
                             '63 - IMEXRK2(2,3,2), 64 - IMEX_DIRK_RK3.')
    parser.add_argument('-tf', '--t-final', default=10.0, type=float,
                        help='Final time; start time is 0.')
    parser.add_argument('-dt', '--time-step', default=0.01, type=float,
                        help='Time step.')
    parser.add_argument('-dc', '--diffusion-coeff', default=0.01, type=float,
                        help='Diffusion coefficient in the PDE.')
    parser.add_argument('-cg', '--continuous-galerkin', action='store_true',
                        help='Use Continuous-Galerkin finite elements (default is DG).')
    parser.add_argument('-vs', '--visualization-steps', default=50, type=int,
                        help='Visualize every n-th timestep.')
    parser.add_argument('-no-vis', '--no-visualization', action='store_true',
                        help='Disable GLVis visualization.')
    parser.add_argument('-visit', '--visit-datafiles', action='store_true',
                        help='Save ASCII data files for VisIt visualization.')
    parser.add_argument('-paraview', '--paraview-datafiles', action='store_true',
                        help='Save data files for ParaView visualization.')
    args = parser.parse_args()
    parser.print_options(args)
    meshfile = expanduser(join(os.path.dirname(__file__), '..', 'data', args.mesh))
    run(meshfile=meshfile, problem=args.problem, ref_levels=args.refine,
        order=args.order, ode_solver_type=args.ode_solver, t_final=args.t_final,
        dt=args.time_step, diffusion_term=args.diffusion_coeff,
        cg=args.continuous_galerkin, vis_steps=args.visualization_steps,
        visualization=not args.no_visualization, visit=args.visit_datafiles,
        paraview=args.paraview_datafiles)
