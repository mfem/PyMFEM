"""Standalone serial DG advection example with native NSB callbacks.

All coefficients and the explicit evolution operator use public MFEM directors.
Supports CPU full assembly and the explicit ODE schemes of examples/ex9.py.
Scratch state is reusable: callbacks on each instance must not be concurrent
or reentrant. Matrix, solver and backing-array owners stay in Python setup.
"""

import sys
import math
from pathlib import Path
import numpy as np
from numba import njit, types
import nsb_rt as nsb
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def initial_value(x, bb_min, bb_max, problem):
    dim = x.size
    center = (bb_min + bb_max) / 2.0
    X = 2.0 * (x - center) / (bb_max - bb_min)
    if problem == 0 or problem == 1:
        if dim == 1:
            return math.exp(-40.0 * (X[0] - 0.5)**2)
        rx, ry = 0.45, 0.25
        cx, cy, w = 0.0, -0.2, 10.0
        if dim == 3:
            s = 1.0 + 0.25 * math.cos(2.0 * math.pi * x[2])
            rx *= s
            ry *= s
        return (math.erfc(w * (X[0]-cx-rx)) *
                math.erfc(-w * (X[0]-cx+rx)) *
                math.erfc(w * (X[1]-cy-ry)) *
                math.erfc(-w * (X[1]-cy+ry))) / 16.0
    if problem == 2:
        rho = math.hypot(x[0], x[1])
        phi = math.atan2(x[1], x[0])
        return math.sin(math.pi * rho)**2 * math.sin(3.0 * phi)
    return math.sin(math.pi * X[0]) * math.sin(math.pi * X[1])


def get_directors():
    """Compile the directors after checking the execution mode and mesh."""
    COEFFICIENT_STATE = {
        "point": mfem.Vector, "bb_min": types.float64[::1],
        "bb_max": types.float64[::1], "problem": types.intc,
    }


    @nsb.director(state=COEFFICIENT_STATE, fallback="silent")
    class velocity_coeff(mfem.VectorCoefficient):
        def __init__(self, bb_min, bb_max, problem):
            super().__init__(len(bb_min))
            self.point = mfem.Vector(len(bb_min))
            self.bb_min = np.ascontiguousarray(bb_min, dtype=np.float64)
            self.bb_max = np.ascontiguousarray(bb_max, dtype=np.float64)
            self.problem = problem

        @nsb.override(stateaccess="*")
        def Eval(self, V, T, ip):
            T.Transform(ip, self.point)
            x = self.point.GetDataArray()
            dim = x.size
            V.SetSize(dim)
            v = V.GetDataArray()
            if dim == 1:
                v[0] = 1.0
            elif self.problem == 0:
                if dim == 2:
                    v[0], v[1] = math.sqrt(2.0/3.0), math.sqrt(1.0/3.0)
                else:
                    v[0] = math.sqrt(3.0/6.0)
                    v[1] = math.sqrt(2.0/6.0)
                    v[2] = math.sqrt(1.0/6.0)
            else:
                center = (self.bb_min + self.bb_max) / 2.0
                X = 2.0 * (x - center) / (self.bb_max - self.bb_min)
                w = math.pi / 2.0
                if self.problem == 3:
                    d = max((X[0]+1.0)*(1.0-X[0]), 0.0) * max((X[1]+1.0)*(1.0-X[1]), 0.0)
                    w *= d*d
                v[0], v[1] = w*X[1], -w*X[0]
                if dim == 3:
                    v[2] = 0.0


    @nsb.director(state=COEFFICIENT_STATE, fallback="silent")
    class u0_coeff(mfem.Coefficient):
        def __init__(self, bb_min, bb_max, problem):
            super().__init__()
            self.point = mfem.Vector(len(bb_min))
            self.bb_min = np.ascontiguousarray(bb_min, dtype=np.float64)
            self.bb_max = np.ascontiguousarray(bb_max, dtype=np.float64)
            self.problem = problem

        @nsb.override(stateaccess="*")
        def Eval(self, T, ip):
            T.Transform(ip, self.point)
            return initial_value(self.point.GetDataArray(), self.bb_min, self.bb_max, self.problem)


    @nsb.director(state={"value": types.float64}, fallback="silent")
    class inflow_coeff(mfem.Coefficient):
        def __init__(self):
            super().__init__()
            self.value = 0.0

        @nsb.override
        def Eval(self, T, ip):
            return self.value


    @nsb.director(state={"M": mfem.SparseMatrix, "K": mfem.SparseMatrix,
                         "b": mfem.Vector, "z": mfem.Vector,
                         "M_solver": mfem.CGSolver}, fallback="silent")
    class FE_Evolution(mfem.TimeDependentOperator):
        def __init__(self, M, K, b, *, owners=()):
            super().__init__(M.Size())
            self.M = M
            self.K = K
            # LinearForm inherits Vector; use an explicit ordinary Vector proxy.
            self.b = mfem.Vector(b.GetDataArray())
            self.z = mfem.Vector(M.Size())
            self.M_prec = mfem.DSmoother()
            self.M_solver = mfem.CGSolver()
            self.M_solver.SetPreconditioner(self.M_prec)
            self.M_solver.SetOperator(M)
            self.M_solver.iterative_mode = False
            self.M_solver.SetRelTol(1e-9)
            self.M_solver.SetAbsTol(0.0)
            self.M_solver.SetMaxIter(100)
            self.M_solver.SetPrintLevel(0)
            # SpMat() and solver/preconditioner connections borrow their owners.
            self._owners = (*owners, b, self.M_prec)

        @nsb.override(stateaccess="*")
        def Mult(self, x, y):
            self.K.Mult(x, self.z)
            z = self.z.GetDataArray()
            b = self.b.GetDataArray()
            for index in range(z.size):
                z[index] += b[index]
            self.M_solver.Mult(self.z, y)


    return velocity_coeff, u0_coeff, inflow_coeff, FE_Evolution


def run(ref_levels=2, order=3, ode_solver_type=4, t_final=10.0, dt=0.01,
        problem=0, vis_steps=5, meshfile="", visualization=False,
        device="cpu", pa=False, paraview=False):

    if pa or device != "cpu":
        raise ValueError("bridge ex9 supports CPU full assembly only")

    if problem not in (0, 1, 2, 3) or dt <= 0 or t_final < 0 or vis_steps <= 0:
        raise ValueError("invalid problem, time interval or visualization frequency")

    device = mfem.Device(device)
    device.Print()

    # 2. Read the mesh from the given mesh file. We can handle geometrically
    #    periodic meshes in this code.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    if mesh.SpaceDimension() != dim or dim not in (1, 2, 3):
        raise ValueError("bridge ex9 requires a full-dimensional mesh in 1D, 2D or 3D")
    if dim == 1 and problem in (2, 3):
        raise ValueError("problems 2 and 3 require dimension at least two")

    # 3. Define the ODE solver used for time integration. Several explicit
    #    Runge-Kutta methods are available.
    factories = {1: mfem.ForwardEulerSolver, 2: lambda: mfem.RK2Solver(1.0),
                 3: mfem.RK3SSPSolver, 4: mfem.RK4Solver, 6: mfem.RK6Solver}
    if ode_solver_type not in factories:
        print("Unknown ODE solver type: " + str(ode_solver_type))

        sys.exit()

    ode_solver = factories[ode_solver_type]()

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement, where 'ref_levels' is a
    #    command-line parameter. If the mesh is of NURBS type, we convert it to
    #    a (piecewise-polynomial) high-order mesh.
    for lev in range(ref_levels):
        mesh.UniformRefinement()
        if mesh.NURBSext:
            mesh.SetCurvature(max(order, 1))
    bb_min, bb_max = mesh.GetBoundingBox(max(order, 1))

    if np.any(np.asarray(bb_max) <= np.asarray(bb_min)):
        raise ValueError("mesh bounding box must have positive extent")

    # 5. Define the discontinuous DG finite element space of the given
    #    polynomial order on the refined mesh.
    fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
    fes = mfem.FiniteElementSpace(mesh, fec)
    print("Number of unknowns:", fes.GetVSize())

    # 6. Set up and assemble the bilinear and linear forms corresponding to the
    #    DG discretization. The DGTraceIntegrator involves integrals over mesh
    #    interior faces.
    velocity_coeff, u0_coeff, inflow_coeff, FE_Evolution = get_directors()

    velocity = velocity_coeff(bb_min, bb_max, problem)
    inflow = inflow_coeff()
    u0 = u0_coeff(bb_min, bb_max, problem)
    adv = None

    # The try/finally block ensures proper director cleanup. See the note in
    # the finally section for details.
    try:
        # 6. Assemble the bilinear and linear forms for the DG discretization.
        m = mfem.BilinearForm(fes)
        m.AddDomainIntegrator(mfem.MassIntegrator())

        k = mfem.BilinearForm(fes)
        k.AddDomainIntegrator(mfem.ConvectionIntegrator(velocity, -1.0))
        k.AddInteriorFaceIntegrator(
            mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5)))
        k.AddBdrFaceIntegrator(
            mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5)))

        b = mfem.LinearForm(fes)
        b.AddBdrFaceIntegrator(mfem.BoundaryFlowIntegrator(inflow, velocity, -1.0, -0.5))

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
        mesh.Print("ex9.mesh", 8)
        u.Save("ex9-init.gf", 8)

        if paraview:
            pd = mfem.ParaViewDataCollection("Example9", mesh)
            pd.SetPrefixPath("ParaView")
            pd.RegisterField("solution", u)
            pd.SetLevelsOfDetail(order)
            pd.SetDataFormat(mfem.VTKFormat_BINARY)
            pd.SetHighOrderOutput(True)
            pd.SetCycle(0)
            pd.SetTime(0.0)
            pd.Save()

        if visualization:
            sol_sock = mfem.socketstream("localhost", 19916)
            sol_sock.precision(8)
            sol_sock.send_solution(mesh, u)

        adv = FE_Evolution(m.SpMat(), k.SpMat(), b, owners=(m, k))

        ode_solver.Init(adv)
        t, ti = 0.0, 0

        while t <= t_final - dt/2:
            t, dt = ode_solver.Step(u, t, dt)
            ti += 1

            if ti % vis_steps == 0:
                print("time step:", ti, ", time:", np.round(t, 3))
                if paraview:
                    pd.SetCycle(ti)
                    pd.SetTime(t)
                    pd.Save()
                if visualization:
                    sol_sock.send_solution(mesh, u)

        u.Save("ex9-final.gf", 8)
        return {"steps": ti, "time": t, "solution": u.GetDataArray().copy()}

    finally:
        # Call close() to detach the directors and release their callback/state
        # resources.

        # Note:
        #  1) Demonstrate explicit director cleanup, including on exceptions.
        #  2) The cleanup order follows these lifetime requirements:
        #    * Release the ODE solver before closing adv, because it holds an operator
        #      reference.
        #    * Keep matrices and the mass preconditioner alive while adv can run.
        #    * Keep coefficient directors alive while forms may assemble or grid functions
        #      may project them.
        #  3) On normal completion, these close() calls can be omitted if no
        #     further native use occurs and the Python owners remain alive until
        #     native use ends. Cleanup then occurs when the objects are finalized;
        #     explicit close() releases the resources at a predictable point.

        ode_solver = None
        if adv is not None:
            adv.close()
        velocity.close()
        inflow.close()
        u0.close()

if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser
    parser = ArgParser(description="Ex9 with native NSB callbacks (CPU full assembly)")
    parser.add_argument('-m', '--mesh',
                        default='periodic-hexagon.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-r", "--refine",
                        default=2,
                        action='store', type=int,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument("-p", "--problem",
                        default=0, type=int,
                        help="Problem setup to use. See options in velocity_function().")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-o', '--order',
                        action='store', default=3, type=int,
                        help="Order (degree) of the finite elements.")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    parser.add_argument("-s", "--ode-solver",
                        default=4, type=int,
                        help="\n".join(["ODE solver: 1 - Forward Euler,\n\t",
                                        "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6,\n\t",
                                        "            11 - Backward Euler,\n\t",
                                        "            12 - SDIRK23 (L-stable), 13 - SDIRK33,\n\t",
                                        "            22 - Implicit Midpoint Method,\n\t",
                                        "            23 - SDIRK23 (A-stable), 24 - SDIRK34", ]))
    parser.add_argument("-tf", "--t-final",
                        default=10, type=float,
                        help="Final time; start time is 0.")
    parser.add_argument("-dt", "--time-step",
                        default=0.01, type=float,
                        help="Time step.")
    parser.add_argument("-vs", "--visualization-steps",
                        default=5, type=int,
                        help="Visualize every n-th timestep.")
    parser.add_argument("-paraview", "--paraview-datafiles",
                        default=False,
                        action='store_true',
                        help="Save data files for ParaView (paraview.org) visualization.")

    args = parser.parse_args()
    parser.print_options(args)
    order = args.order
    meshfile = Path(args.mesh).expanduser()
    if not meshfile.is_absolute():
        meshfile = Path(__file__).resolve().parents[2] / "data" / meshfile
    visualization = args.visualization
    device = args.device
    pa = args.partial_assembly
    ode_solver_type = args.ode_solver
    t_final = args.t_final
    dt = args.time_step
    vis_step = args.visualization_steps
    ref_levels = args.refine
    problem = args.problem
    paraview = args.paraview_datafiles

    run(ref_levels=ref_levels,
        order=order,
        ode_solver_type=ode_solver_type,
        t_final=t_final,
        dt=dt,
        problem=problem,
        vis_steps=vis_step,
        meshfile=str(meshfile),
        visualization=visualization,
        device=device,
        pa=pa,
        paraview=paraview)
