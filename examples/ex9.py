'''
   MFEM example 9
      This is a version of Example 1 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
import os
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname, exists
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc


# Inflow boundary condition (zero for the problems considered in this example)
class inflow_coeff(mfem.PyCoefficient):
    def EvalValue(self, x):
        return 0


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
        self.K.Mult(x, self.z)
        self.z += self.b
        self.M_solver.Mult(self.z, y)


def run(ref_levels=2,
        order=3,
        ode_solver_type=4,
        t_final=10,
        dt=0.01,
        problem=0,
        vis_steps=5,
        meshfile='',
        visualization=False,
        device='cpu',
        pa=False,
        paraview=False):

    device = mfem.Device(device)
    device.Print()
    
    # 2. Read the mesh from the given mesh file. We can handle geometrically
    #    periodic meshes in this code.

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    # 3. Define the ODE solver used for time integration. Several explicit
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

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement, where 'ref_levels' is a
    #    command-line parameter. If the mesh is of NURBS type, we convert it to
    #    a (piecewise-polynomial) high-order mesh.
    for lev in range(ref_levels):
        mesh.UniformRefinement()
        if mesh.NURBSext:
            mesh.SetCurvature(max(order, 1))
    bb_min, bb_max = mesh.GetBoundingBox(max(order, 1))

    # 5. Define the discontinuous DG finite element space of the given
    #    polynomial order on the refined mesh.
    fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
    fes = mfem.FiniteElementSpace(mesh, fec)

    print("Number of unknowns: " + str(fes.GetVSize()))

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
                            erfc(w * (X[1]-cy-ry)) * erfc(-w*(X[1]-cy+ry)))/16.

            elif problem == 2:
                rho = hypot(x[0], x[1])
                phi = arctan2(x[1], x[0])
                return (sin(pi * rho) ** 2) * sin(3*phi)
            elif problem == 3:
                return sin(pi * X[0]) * sin(pi * X[1])

            return 0.0

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

    mesh.Print('ex9.mesh', 8)
    u.Save('ex9-init.gf', 8)

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
        sol_sock.send_solution(mesh,  u)

    adv = FE_Evolution(m.SpMat(), k.SpMat(), b)

    ode_solver.Init(adv)
    t = 0.0
    ti = 0
    while True:
        if t > t_final - dt/2:
            break
        t, dt = ode_solver.Step(u, t, dt)
        ti = ti + 1

        if ti % vis_steps == 0:
            print("time step: " + str(ti) + ", time: " + str(np.round(t, 3)))
            if paraview:
                pd.SetCycle(ti)
                pd.SetTime(t)
                pd.Save()
            if visualization:
                sol_sock.send_solution(mesh,  u)

    u.Save('ex9-final.gf', 8)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex9')
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
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
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
        meshfile=meshfile,
        visualization=visualization,
        device=device,
        pa=pa,
        paraview=paraview)
