'''
   MFEM example 18
      This is a version of Example 18 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
from ex18_common import FE_Evolution, InitialCondition, RiemannSolver, DomainIntegrator, FaceIntegrator
from mfem.common.arg_parser import ArgParser
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

# Equation constant parameters.(using globals to share them with ex18_common)
import ex18_common


def run(problem=1,
        ref_levels=1,
        order=3,
        ode_solver_type=4,
        t_final=0.5,
        dt=-0.01,
        cfl=0.3,
        visualization=True,
        vis_steps=50,
        meshfile=''):

    ex18_common.num_equation = 4
    ex18_common.specific_heat_ratio = 1.4
    ex18_common.gas_constant = 1.0
    ex18_common.problem = problem
    num_equation = ex18_common.num_equation

    # 2. Read the mesh from the given mesh file. This example requires a 2D
    #    periodic mesh, such as ../data/periodic-square.mesh.
    meshfile = expanduser(join(dirname(__file__), '..', 'data', meshfile))
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
    #    command-line parameter.
    for lev in range(ref_levels):
        mesh.UniformRefinement()

    # 5. Define the discontinuous DG finite element space of the given
    #    polynomial order on the refined mesh.

    fec = mfem.DG_FECollection(order, dim)
    # Finite element space for a scalar (thermodynamic quantity)
    fes = mfem.FiniteElementSpace(mesh, fec)
    # Finite element space for a mesh-dim vector quantity (momentum)
    dfes = mfem.FiniteElementSpace(mesh, fec, dim, mfem.Ordering.byNODES)
    # Finite element space for all variables together (total thermodynamic state)
    vfes = mfem.FiniteElementSpace(
        mesh, fec, num_equation, mfem.Ordering.byNODES)

    assert fes.GetOrdering() == mfem.Ordering.byNODES, "Ordering must be byNODES"
    print("Number of unknowns: " + str(vfes.GetVSize()))

    # 6. Define the initial conditions, save the corresponding mesh and grid
    #    functions to a file. This can be opened with GLVis with the -gc option.
    #    The solution u has components {density, x-momentum, y-momentum, energy}.
    #    These are stored contiguously in the BlockVector u_block.

    offsets = [k*vfes.GetNDofs() for k in range(num_equation+1)]
    offsets = mfem.intArray(offsets)
    u_block = mfem.BlockVector(offsets)
    mom = mfem.GridFunction(dfes, u_block,  offsets[1])

    #
    #  Define coefficient using VecotrPyCoefficient and PyCoefficient
    #  A user needs to define EvalValue method
    #
    u0 = InitialCondition(num_equation)
    sol = mfem.GridFunction(vfes, u_block.GetData())
    sol.ProjectCoefficient(u0)

    mesh.Print("vortex.mesh", 8)
    for k in range(num_equation):
        uk = mfem.GridFunction(fes, u_block.GetBlock(k).GetData())
        sol_name = "vortex-" + str(k) + "-init.gf"
        uk.Save(sol_name, 8)

    #  7. Set up the nonlinear form corresponding to the DG discretization of the
    #     flux divergence, and assemble the corresponding mass matrix.
    Aflux = mfem.MixedBilinearForm(dfes, fes)
    Aflux.AddDomainIntegrator(DomainIntegrator(dim))
    Aflux.Assemble()

    A = mfem.NonlinearForm(vfes)
    rsolver = RiemannSolver()
    ii = FaceIntegrator(rsolver, dim)
    A.AddInteriorFaceIntegrator(ii)

    #  8. Define the time-dependent evolution operator describing the ODE
    #     right-hand side, and perform time-integration (looping over the time
    #     iterations, ti, with a time-step dt).
    euler = FE_Evolution(vfes, A, Aflux.SpMat())

    if (visualization):
        sout = mfem.socketstream("localhost", 19916)
        sout.precision(8)
        sout << "solution\n" << mesh << mom
        sout << "pause\n"
        sout.flush()
        print("GLVis visualization paused.")
        print(" Press space (in the GLVis window) to resume it.")

    # Determine the minimum element size.
    hmin = 0
    if (cfl > 0):
        hmin = min([mesh.GetElementSize(i, 1) for i in range(mesh.GetNE())])

    t = 0.0
    euler.SetTime(t)
    ode_solver.Init(euler)
    if (cfl > 0):
        #  Find a safe dt, using a temporary vector. Calling Mult() computes the
        #  maximum char speed at all quadrature points on all faces.
        z = mfem.Vector(A.Width())
        A.Mult(sol, z)

        dt = cfl * hmin / ex18_common.max_char_speed / (2*order+1)

    # Integrate in time.
    done = False
    ti = 0
    while not done:
        dt_real = min(dt, t_final - t)

        t, dt_real = ode_solver.Step(sol, t, dt_real)

        if (cfl > 0):
            dt = cfl * hmin / ex18_common.max_char_speed / (2*order+1)
        ti = ti+1
        done = (t >= t_final - 1e-8*dt)
        if (done or ti % vis_steps == 0):
            print("time step: " + str(ti) + ", time: " + "{:g}".format(t))
            if (visualization):
                sout << "solution\n" << mesh << mom << flush

    #  9. Save the final solution. This output can be viewed later using GLVis:
    #    "glvis -m vortex.mesh -g vortex-1-final.gf".
    for k in range(num_equation):
        uk = mfem.GridFunction(fes, u_block.GetBlock(k).GetData())
        sol_name = "vortex-" + str(k) + "-final.gf"
        uk.Save(sol_name, 8)

    print(" done")
    # 10. Compute the L2 solution error summed for all components.
    if (t_final == 2.0):
        error = sol.ComputeLpError(2., u0)
        print("Solution error: " + str(error))


if __name__ == "__main__":

    parser = ArgParser(description='Ex18')
    parser.add_argument('-m', '--mesh',
                        default='periodic-square.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-p', '--problem',
                        action='store', default=1, type=int,
                        help='Problem setup to use. See options in velocity_function().')
    parser.add_argument('-r', '--refine',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly.")
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

    parser.print_options(args)

    run(problem=args.problem,
        ref_levels=args.refine,
        order=args.order,
        ode_solver_type=args.ode_solver,
        t_final=args.t_final,
        dt=args.time_step,
        cfl=args.cfl_number,
        visualization=args.visualization,
        vis_steps=args.visualization_steps,
        meshfile=args.mesh)
