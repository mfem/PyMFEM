#  MFEM Example 41 (converted from ex40.cpp)
# 
#  Sample runs:
#   python ex41.py
#   python ex41.py -cg
#   python ex41.py -m ../data/periodic-hexagon.mesh -p 0 -r 2 -dt 0.005 -tf 10
#   python ex41.py -m ../data/periodic-square.mesh -p 1 -r 2 -dt 0.005 -tf 9
#   python ex41.py -m ../data/periodic-hexagon.mesh -p 1 -r 2 -dt 0.005 -tf 9
#   python ex41.py -m ../data/amr-quad.mesh -p 1 -r 2 -dt 0.002 -tf 9
#   python ex41.py -m ../data/star-q3.mesh -p 1 -r 2 -dt 0.001 -tf 9
#   python ex41.py -m ../data/star-mixed.mesh -p 1 -r 2 -dt 0.005 -tf 9
#   python ex41.py -m ../data/disc-nurbs.mesh -p 1 -r 3 -dt 0.005 -tf 9
#   python ex41.py -m ../data/disc-nurbs.mesh -p 2 -r 3 -dt 0.005 -tf 9
#   python ex41.py -m ../data/periodic-square.mesh -p 3 -r 4 -dt 0.0025 -tf 9 -vs 20
#   python ex41.py -m ../data/periodic-cube.mesh -p 0 -r 2 -o 2 -dt 0.01 -tf 8
# 
#  Device sample runs:
# 
#  Description:  This example code solves the time-dependent advection-diffusion
#                equation du/dt + v.grad(u) - a div(grad(u)) = 0, where v is a
#                given fluid velocity, a is the diffusion coefficient, and
#                u0(x)=u(0,x) is a given initial condition.
# 
#                The example demonstrates the use of Discontinuous Galerkin (DG)
#                bilinear forms in MFEM (face integrators), and the use of IMEX
#                ODE time integrators.
# 
#                The option to use continuous finite elements is available too.


from numpy import zeros
from math import erfc, pi, sin, cos, erfc, atan2, exp
from sys.float_info import epsilon

from mfem import ser as mfem
from mfem.common.arg_parser import ArgParser


# Global bounding box coordinates.
bb_min = None
bb_max = None

# Velocity coefficient.
def make_velocity(problem, dim, bb_min, bb_max):

    if dim == 1:
        @mfem.jit.vector(sdim=1, vdim=1)
        def velocity(x):
            v = zeros(1)
            v[0] = 1.0
            return v

    elif dim == 2:

        @mfem.jit.vector(sdim=2, vdim=2)
        def velocity(x):

            X = zeros(2)

            for i in range(2):
                center = 0.5 * (bb_min[i] + bb_max[i])
                X[i] = 2.0 * (x[i] - center) / \
                       (bb_max[i] - bb_min[i])

            v = zeros(2)

            if problem == 0:
                # translation
                v[0] = sqrt(2.0 / 3.0)
                v[1] = sqrt(1.0 / 3.0)

            elif problem in [1, 2]:
                # rigid rotation
                w = pi / 2.0

                v[0] = w * X[1]
                v[1] = -w * X[0]

            elif problem == 3:
                # Clockwise twisting rotation in 2D around the origin
                w = pi / 2.0

                d = max((X[0] + 1.0) * (1.0 - X[0]), 0.0)
                d *= max((X[1] + 1.0) * (1.0 - X[1]), 0.0)
                d = d * d

                v[0] = d * w * X[1]
                v[1] = -d * w * X[0]

            return v

    else:
        @mfem.jit.vector(sdim=3, vdim=3)
        def velocity(x):

            X = zeros(3)

            for i in range(3):
                center = 0.5 * (bb_min[i] + bb_max[i])
                X[i] = 2.0 * (x[i] - center) / \
                       (bb_max[i] - bb_min[i])

            v = zeros(3)

            if problem == 0:
                v[0] = sqrt(3.0 / 6.0)
                v[1] = sqrt(2.0 / 6.0)
                v[2] = sqrt(1.0 / 6.0)
                
            elif problem in [1, 2]:
                w = pi / 2.0
                
                v[0] = w * X[1]
                v[1] = - w * X[0]
                v[2] = 0.0
                
            else:
                # Clockwise twisting rotation in 2D around the origin                
                w = pi / 2.0

                d = max((X[0] + 1.0) * (1.0 - X[0]), 0.0)
                d *= max((X[1] + 1.0) * (1.0 - X[1]), 0.0)
                d = d * d

                v[0] = d * w * X[1]
                v[1] = -d * w * X[0]
                v[2] = 0.0

            return v

    return velocity

# Initial condition coefficient.
def make_initial_condition(problem, dim, bb_min, bb_max):

    @mfem.jit.scalar(sdim=dim)
    def u0(x):

        X = zeros(dim)

        for i in range(dim):
            center = 0.5 * (bb_min[i] + bb_max[i])
            X[i] = 2.0 * (x[i] - center) / (bb_max[i] - bb_min[i])

        if problem in [0, 1]:
            if dim == 1:
                return exp(-40.0 * (X[0] - 0.5)**2)

            rx = 0.45
            ry = 0.25
            cx = 0.0
            cy = -0.2
            w = 10.

            if dim == 3:
                rx *= (1. + 0.25 * cos( 2 * pi * X[2]))
                ry *= (1. + 0.25 * cos( 2 * pi * X[2]))
                
            return (erfc(w*(X[0]-cx-rx))*erfc(-w*(X[0]-cx+rx)) *
                    erfc(w*(X[1]-cy-ry))*erfc(-w*(X[1]-cy+ry)) )/16
                
        elif problem == 2:
            rho = sqrt(x[0]**2 + y[2]**2)
            phi = atan2(y[1], x[0])
            return sin(pi * rho)**2 * sin(3 * phi)

        elif problem == 3:
            return sin(pi * X[0]) * sin(pi * X[1])

        return 0.0

    return u0


# -----------------------------------------------------------------------------
# Time-dependent operator
# -----------------------------------------------------------------------------

# Solver for the implicit part of the ODE (the diffusion term).
# Solves systems of the form: (M + dt*S) k = rhs.
class ImplicitSolver(mfem.Solver):

    def __init__(self, M, S, fes):

        mfem.PySolver.__init__(self, M.Height(), M.Width())

        self.M = M
        self.S = S
        self.A = None

        self.prec = mfem.DSmoother()
        self.solver = mfem.CGSolver()
        self.dt = 1.0

        self.solver.iterative_mode = False
        self.solver.SetRelTol(1e-9)
        self.solver.SetAbsTol(0.0)
        self.solver.SetMaxIter(100)
        self.solver.SetPreconditioner(self.prec)

    def SetTimeStep(self, dt_):
        ddt = self.dt - dt_

        if abs(ddt) > epsilon*10:
            # Form operator A = M + dt*S
            self.A = mfem.Add(1.0, self.M, dt, self.S)
            self.solver.SetOperator(self.A)

#   A time-dependent operator for the right-hand side of the ODE. The weak
#   form of the advection-diffusion equation is M du/dt = K u - S u + b,
#   where M is the mass matrix, K and S are the advection and diffusion
#   matrices, and b describes the flow on the boundary. In the case of IMEX
#   evolution, the diffusion term is treated implicitly, and the advection
#   term is treated explicitly.  */
class IMEX_Evolution(mfem.PyTimeDependentOperator):

    def __init__(self,
                 M,
                 K,
                 S):

        mfem.PyTimeDependentOperator.__init__(
            self,
            M.Height()
        )

        self.M = M
        self.K = K
        self.S = S

        self.z = mfem.Vector(M.Height())
        self.rhs = mfem.Vector(M.Height())

        # -------------------------------------------------------------
        # Mass matrix inverse
        # -------------------------------------------------------------

        self.M_prec = mfem.DSmoother()

        self.M_solver = mfem.CGSolver()

        self.M_solver.SetRelTol(1e-9)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(200)
        self.M_solver.SetPrintLevel(0)

        self.M_solver.SetPreconditioner(
            self.M_prec
        )

        self.M_solver.SetOperator(M)

        # -------------------------------------------------------------
        # Implicit diffusion solver
        # -------------------------------------------------------------

        self.implicit_solver = ImplicitSolver(
            M,
            S
        )

    # -----------------------------------------------------------------
    # Explicit part
    #
    #     y = M^{-1} (-K x)
    # -----------------------------------------------------------------

    def Mult(self, x, y):

        self.K.Mult(x, self.z)

        self.z *= -1.0

        self.M_solver.Mult(self.z, y)

    # -----------------------------------------------------------------
    # Implicit solve
    #
    # Solve:
    #
    #     (M + dt S) k = -S x
    # -----------------------------------------------------------------

    def ImplicitSolve(self,
                      dt,
                      x,
                      k):

        self.S.Mult(x, self.rhs)

        self.rhs *= -1.0

        self.implicit_solver.SetTimeStep(dt)

        self.implicit_solver.Mult(
            self.rhs,
            k
        )


# -----------------------------------------------------------------------------
# Time-dependent evolution operator.
#
# This class implements the semi-discrete operator:
#
#     M du/dt = -(K + D) u
#
# where:
#
#     M : mass matrix
#     K : advection operator
#     D : diffusion operator
#
# The operator supports both explicit and implicit time integration.
# -----------------------------------------------------------------------------


class FE_Evolution(mfem.PyTimeDependentOperator):

    def __init__(self, M, K, D, dt, implicit):

        mfem.PyTimeDependentOperator.__init__(self, M.Height())

        self.M = M
        self.K = K
        self.D = D
        self.dt = dt
        self.implicit = implicit

        self.z = mfem.Vector(M.Height())
        self.tmp = mfem.Vector(M.Height())

        # Mass solver
        self.M_prec = mfem.DSmoother()

        self.M_solver = mfem.CGSolver()
        self.M_solver.SetRelTol(1e-9)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(200)
        self.M_solver.SetPrintLevel(0)

        self.M_solver.SetPreconditioner(
            self.M_prec
        )

        self.M_solver.SetOperator(M)

        # Implicit operator
        if implicit:

            self.T = mfem.Add(1.0, M, dt, D)

            self.T_prec = mfem.DSmoother()

            self.T_solver = mfem.CGSolver()
            self.T_solver.SetRelTol(1e-9)
            self.T_solver.SetAbsTol(0.0)
            self.T_solver.SetMaxIter(200)
            self.T_solver.SetPrintLevel(0)

            self.T_solver.SetPreconditioner(
                self.T_prec
            )

            self.T_solver.SetOperator(self.T)

    # -----------------------------------------------------------------

    def Mult(self, x, y):

        self.K.Mult(x, self.z)

        self.D.Mult(x, self.tmp)

        self.z += self.tmp

        self.z *= -1.0

        self.M_solver.Mult(self.z, y)

    # -----------------------------------------------------------------

    def ImplicitSolve(self, dt, x, k):

        self.K.Mult(x, self.z)

        self.z *= -1.0

        self.T_solver.Mult(self.z, k)


# -----------------------------------------------------------------------------
# Main program.
#
# Parse command-line arguments and launch the simulation.
# -----------------------------------------------------------------------------


# Main driver routine.
#
# This routine:
#
#   1. Creates and refines the mesh
#   2. Constructs the finite element space
#   3. Builds the advection and diffusion operators
#   4. Initializes the ODE solver
#   5. Evolves the solution in time
#   6. Sends output to GLVis / ParaView / VisIt
# -----------------------------------------------------------------------------


# Main driver
# -----------------------------------------------------------------------------


def run(order,
        ref_levels,
        ode_solver_type,
        t_final,
        dt,
        problem,
        mesh_file,
        visualization,
        vis_steps,
        diffusion,
        cg,
        pa,
        device_config,
        paraview,
        visit):

    global bb_min
    global bb_max

    # -----------------------------------------------------------------
    # Device
    # -----------------------------------------------------------------

    device = mfem.Device(device_config)
    device.Print()

    # -----------------------------------------------------------------
    # Read the mesh from the given mesh file.
    #
    # The mesh can be periodic and can contain high-order curved
    # elements.
    # -----------------------------------------------------------------

    # Mesh
    # -----------------------------------------------------------------

    mesh = mfem.Mesh(mesh_file, 1, 1)

    dim = mesh.Dimension()

    for _ in range(ref_levels):
        mesh.UniformRefinement()

    if mesh.NURBSext:
        mesh.SetCurvature(max(order, 1))

    bb_min_v = mfem.Vector(dim)
    bb_max_v = mfem.Vector(dim)

    mesh.GetBoundingBox(bb_min_v, bb_max_v, max(order, 1))

    bb_min = bb_min_v.GetDataArray()
    bb_max = bb_max_v.GetDataArray()

    # -----------------------------------------------------------------
    # Define the finite element discretization.
    #
    # Continuous H1 elements are used when cg=True.
    # Otherwise discontinuous Galerkin elements are used.
    # -----------------------------------------------------------------

    # Finite element space
    # -----------------------------------------------------------------

    if cg:
        fec = mfem.H1_FECollection(order, dim)
    else:
        fec = mfem.DG_FECollection(
            order,
            dim,
            mfem.BasisType.GaussLobatto
        )

    fes = mfem.FiniteElementSpace(mesh, fec)

    print(f"Number of unknowns: {fes.GetVSize()}")

    # -----------------------------------------------------------------
    # Coefficients
    # -----------------------------------------------------------------

    velocity = make_velocity(
        problem,
        dim,
        bb_min,
        bb_max
    )

    u0 = make_initial_condition(
        problem,
        dim,
        bb_min,
        bb_max
    )

    # -----------------------------------------------------------------
    # Initial condition
    # -----------------------------------------------------------------

    u = mfem.GridFunction(fes)

    u.ProjectCoefficient(u0)

    # -----------------------------------------------------------------
    # Assemble the mass matrix.
    # -----------------------------------------------------------------

    # Mass matrix
    # -----------------------------------------------------------------

    m = mfem.BilinearForm(fes)

    if pa:
        m.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

    m.AddDomainIntegrator(
        mfem.MassIntegrator()
    )

    m.Assemble()
    m.Finalize()

    M = m.SpMat()

    # -----------------------------------------------------------------
    # Assemble the advection operator.
    #
    # In DG mode, numerical fluxes are added through the
    # DGTraceIntegrator.
    # -----------------------------------------------------------------

    # Advection operator
    # -----------------------------------------------------------------

    alpha = -1.0

    k = mfem.BilinearForm(fes)

    if pa:
        k.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

    k.AddDomainIntegrator(
        mfem.ConvectionIntegrator(
            velocity,
            alpha
        )
    )

    if not cg:

        k.AddInteriorFaceIntegrator(
            mfem.TransposeIntegrator(
                mfem.DGTraceIntegrator(
                    velocity,
                    1.0,
                    -0.5
                )
            )
        )

        k.AddBdrFaceIntegrator(
            mfem.TransposeIntegrator(
                mfem.DGTraceIntegrator(
                    velocity,
                    1.0,
                    -0.5
                )
            )
        )

    k.Assemble()
    k.Finalize()

    K = k.SpMat()

    # -----------------------------------------------------------------
    # Assemble the diffusion operator.
    #
    # In DG mode, interior penalty diffusion terms are added.
    # -----------------------------------------------------------------

    # Diffusion operator
    # -----------------------------------------------------------------

    d = mfem.BilinearForm(fes)

    if pa:
        d.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

    diffcoeff = mfem.ConstantCoefficient(diffusion)

    d.AddDomainIntegrator(
        mfem.DiffusionIntegrator(diffcoeff)
    )

    if not cg:

        sigma = -1.0
        kappa = (order + 1) ** 2

        d.AddInteriorFaceIntegrator(
            mfem.DGDiffusionIntegrator(
                diffcoeff,
                sigma,
                kappa
            )
        )

        d.AddBdrFaceIntegrator(
            mfem.DGDiffusionIntegrator(
                diffcoeff,
                sigma,
                kappa
            )
        )

    d.Assemble()
    d.Finalize()

    D = d.SpMat()

    # -----------------------------------------------------------------
    # Evolution operator
    # -----------------------------------------------------------------

    implicit = ode_solver_type in [11, 12, 13, 22, 23]

    adv = FE_Evolution(
        M,
        K,
        D,
        dt,
        implicit
    )

    # -----------------------------------------------------------------
    # Select the ODE time integrator.
    #
    # The numbering follows the MFEM example convention.
    # -----------------------------------------------------------------

    # ODE solver
    # -----------------------------------------------------------------

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

    elif ode_solver_type == 11:
        ode_solver = mfem.BackwardEulerSolver()

    elif ode_solver_type == 12:
        ode_solver = mfem.SDIRK23Solver(2)

    elif ode_solver_type == 13:
        ode_solver = mfem.SDIRK33Solver(3)

    elif ode_solver_type == 22:
        ode_solver = mfem.ImplicitMidpointSolver()

    elif ode_solver_type == 23:
        ode_solver = mfem.SDIRK23Solver(3)

    else:
        raise ValueError("Unknown ODE solver type")

    ode_solver.Init(adv)

    # -----------------------------------------------------------------
    # Visualization
    # -----------------------------------------------------------------

    if visualization:

        sout = mfem.socketstream(
            "localhost",
            19916
        )

        sout.precision(8)

        sout << "solution\n" << mesh << u
        sout << "window_title 'PyMFEM python ex41.py'\n"
        sout.flush()

    # -----------------------------------------------------------------
    # ParaView output
    # -----------------------------------------------------------------

    if paraview:

        pdc = mfem.ParaViewDataCollection(
            "Example41",
            mesh
        )

        pdc.SetPrefixPath("ParaView")
        pdc.RegisterField("solution", u)
        pdc.SetLevelsOfDetail(order)
        pdc.SetDataFormat(mfem.VTKFormat_BINARY)
        pdc.SetHighOrderOutput(True)
        pdc.SetCycle(0)
        pdc.SetTime(0.0)
        pdc.Save()

    # -----------------------------------------------------------------
    # VisIt output
    # -----------------------------------------------------------------

    if visit:

        visit_dc = mfem.VisItDataCollection(
            "Example41",
            mesh
        )

        visit_dc.RegisterField(
            "solution",
            u
        )

        visit_dc.SetCycle(0)
        visit_dc.SetTime(0.0)
        visit_dc.Save()

    # -----------------------------------------------------------------
    # Time integration loop.
    #
    # Advance the solution until the final time is reached.
    # -----------------------------------------------------------------

    # Time integration loop
    # -----------------------------------------------------------------

    t = 0.0
    ti = 0

    while t < t_final - 1e-12:

        dt_real = min(dt, t_final - t)

        ode_solver.Step(u, t, dt_real)

        ti += 1

        if ti % vis_steps == 0 or \
           t >= t_final - 1e-12:

            print(
                f"step = {ti}, "
                f"t = {t:.6e}"
            )

            if visualization:
                sout << "solution\n" << mesh << u
                sout.flush()

            if paraview:
                pdc.SetCycle(ti)
                pdc.SetTime(t)
                pdc.Save()

            if visit:
                visit_dc.SetCycle(ti)
                visit_dc.SetTime(t)
                visit_dc.Save()

    # -----------------------------------------------------------------
    # Save final state
    # -----------------------------------------------------------------

    mesh.Print("python ex41.py.mesh", 8)

    u.Save("python ex41.py.gf", 8)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


if __name__ == "__main__":

    parser = ArgParser(
        description="PyMFEM translation of MFEM python ex41.py"
    )

    parser.add_argument(
        "-m", "--mesh",
        default="../../data/periodic-square.mesh",
        type=str,
        help="Mesh file"
    )

    parser.add_argument(
        "-o", "--order",
        default=3,
        type=int,
        help="Finite element order"
    )

    parser.add_argument(
        "-r", "--refine",
        default=2,
        type=int,
        help="Refinement levels"
    )

    parser.add_argument(
        "-tf", "--t-final",
        default=10.0,
        type=float,
        help="Final time"
    )

    parser.add_argument(
        "-dt", "--time-step",
        default=0.01,
        type=float,
        help="Time step"
    )

    parser.add_argument(
        "-p", "--problem",
        default=0,
        type=int,
        help="Problem type"
    )

    parser.add_argument(
        "-k", "--diffusion",
        default=0.01,
        type=float,
        help="Diffusion coefficient"
    )

    parser.add_argument(
        "-cg", "--continuous",
        action="store_true",
        default=False,
        help="Use H1 space"
    )

    parser.add_argument(
        "-dg", "--discontinuous",
        action="store_false",
        dest="continuous",
        help="Use DG space"
    )

    parser.add_argument(
        "-pa", "--partial-assembly",
        action="store_true",
        default=False,
        help="Enable partial assembly"
    )

    parser.add_argument(
        "-d", "--device",
        default="cpu",
        type=str,
        help="Device configuration"
    )

    parser.add_argument(
        "-s", "--ode-solver",
        default=4,
        type=int,
        help="ODE solver type"
    )

    parser.add_argument(
        "-vis", "--visualization",
        action="store_true",
        default=True,
        help="Enable GLVis visualization"
    )

    parser.add_argument(
        "-no-vis", "--no-visualization",
        action="store_false",
        dest="visualization",
        help="Disable GLVis visualization"
    )

    parser.add_argument(
        "-vs", "--visualization-steps",
        default=50,
        type=int,
        help="Visualization interval"
    )

    parser.add_argument(
        "-visit", "--visit-datafiles",
        action="store_true",
        default=False,
        help="Enable VisIt output"
    )

    parser.add_argument(
        "-paraview", "--paraview-datafiles",
        action="store_true",
        default=False,
        help="Enable ParaView output"
    )

    args = parser.parse_args()

    parser.print_options(args)

    run(order=args.order,
        ref_levels=args.refine,
        ode_solver_type=args.ode_solver,
        t_final=args.t_final,
        dt=args.time_step,
        problem=args.problem,
        mesh_file=args.mesh,
        visualization=args.visualization,
        vis_steps=args.visualization_steps,
        diffusion=args.diffusion,
        cg=args.continuous,
        pa=args.partial_assembly,
        device_config=args.device,
        paraview=args.paraview_datafiles,
        visit=args.visit_datafiles)
