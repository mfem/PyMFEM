'''
   MFEM example 18 parallel
      This is a version of Example 18 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
   Sample runs:

       mpirun -np 4 python ex18p.py -p 1 -rs 2 -rp 1 -o 1 -s 3
       mpirun -np 4 python ex18p.py -p 1 -rs 1 -rp 1 -o 3 -s 4
       mpirun -np 4 python ex18p.py -p 1 -rs 1 -rp 1 -o 5 -s 6
       mpirun -np 4 python ex18p.py -p 2 -rs 1 -rp 1 -o 1 -s 3 -mf
       mpirun -np 4 python ex18p.py -p 2 -rs 1 -rp 1 -o 3 -s 3 -mf

'''
import mfem.par as mfem

from mfem.common.arg_parser import ArgParser

from os.path import expanduser, join, dirname
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc

from ex18_common import (EulerMesh,
                         EulerInitialCondition,
                         DGHyperbolicConservationLaws)

# 1. Initialize MPI.from mpi4py import MPI

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


def run(problem=1,
        ser_ref_levels=0,
        par_ref_levels=1,
        order=3,
        ode_solver_type=4,
        t_final=2.0,
        dt=-0.01,
        cfl=0.3,
        visualization=True,
        vis_steps=50,
        preassembleWeakDiv=False,
        meshfile=''):

    specific_heat_ratio = 1.4
    gas_constant = 1.0
    IntOrderOffset = 1

    device = mfem.Device('cpu')
    if myid == 0:
        device.Print()

    # 2. Read the mesh from the given mesh file. When the user does not provide
    #    mesh file, use the default mesh file for the problem.

    mesh = EulerMesh(meshfile, problem)
    dim = mesh.Dimension()
    num_equation = dim + 2

    # Refine the mesh to increase the resolution. In this example we do
    # 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is a
    # command-line parameter.
    for lev in range(ser_ref_levels):
        mesh.UniformRefinement()

    # Define a parallel mesh by a partitioning of the serial mesh. Refine this
    # mesh further in parallel to increase the resolution. Once the parallel
    # mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh

    # Refine the mesh to increase the resolution. In this example we do
    # 'par_ref_levels' of uniform refinement, where 'par_ref_levels' is a
    # command-line parameter.
    for lev in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 3. Define the ODE solver used for time integration. Several explicit
    #    Runge-Kutta methods are available.
    ode_solver = None
    if ode_solver_type == 1:
        ode_solver = mfem.ForwardEulerSolver()
    elif ode_solver_type == 2:
        ode_solver = mfem.RK2Solver(1.0)
    elif ode_solver_type == 3:
        ode_solver = mfem.RK3SSPSolver()
    elif ode_solver_type == 4:
        ode_solver = mfem.RK4Solver()
    elif ode_solver_type == 6:
        ode_solver = mfem.RK6Solver()
    else:
        print("Unknown ODE solver type: " + str(ode_solver_type))
        exit

    # 4. Define the discontinuous DG finite element space of the given
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

    # 5. Define the initial conditions, save the corresponding mesh and grid
    #    functions to files. These can be opened with GLVis using:
    #    "glvis -np 4  -m euler-mesh -g euler-1-init" (for x-momentum).

    # Initialize the state.
    u0 = EulerInitialCondition(problem,
                               specific_heat_ratio,
                               gas_constant)
    sol = mfem.ParGridFunction(vfes)
    sol.ProjectCoefficient(u0)

    # (Python note): GridFunction pointing to the subset of vector FES.
    #  sol is Vector with dim*fes.GetNDofs()
    #  Since sol.GetDataArray() returns numpy array pointing to the data, we make
    #  Vector from a sub-vector of the returned numpy array and pass it to GridFunction
    #  constructor.
    mom = mfem.GridFunction(dfes, mfem.Vector(
        sol.GetDataArray()[fes.GetNDofs():]))

    # Output the initial solution.
    smyid = '{:0>6d}'.format(myid)
    pmesh.Print("euler-mesh."+smyid, 8)
    for k in range(num_equation):
        uk = mfem.ParGridFunction(fes, mfem.Vector(
            sol.GetDataArray()[k*fes.GetNDofs():]))
        sol_name = "euler-" + str(k) + "-init."+smyid
        uk.Save(sol_name, 8)

    # 6. Set up the nonlinear form with euler flux and numerical flux
    flux = mfem.EulerFlux(dim, specific_heat_ratio)
    numericalFlux = mfem.RusanovFlux(flux)
    formIntegrator = mfem.HyperbolicFormIntegrator(
        numericalFlux, IntOrderOffset)

    euler = DGHyperbolicConservationLaws(vfes, formIntegrator,
                                         preassembleWeakDivergence=preassembleWeakDiv)

    # 7. Visualize momentum with its magnitude
    if (visualization):
        MPI.COMM_WORLD.Barrier()
        sout = mfem.socketstream("localhost", 19916)
        sout.precision(8)
        sout << "parallel " << str(num_procs) << " " << str(myid) << "\n"
        sout << "solution\n" << pmesh << mom
        sout << "window_title 'momentum, t = 0'\n"
        sout << "view 0 0\n"  # view from top
        sout << "keys jlm\n"  # turn off perspective and light, show mesh
        sout << "pause\n"
        sout.flush()

        if myid == 0:
            print("GLVis visualization paused.")
            print(" Press space (in the GLVis window) to resume it.")

    # 8. Time integration
    my_hmin = np.inf
    if (cfl > 0):
        my_hmin = min([pmesh.GetElementSize(i, 1)
                      for i in range(pmesh.GetNE())])

        hmin = MPI.COMM_WORLD.allreduce(my_hmin, op=MPI.MIN)

        # Find a safe dt, using a temporary vector. Calling Mult() computes the
        # maximum char speed at all quadrature points on all faces (and all
        # elements with -mf).
        z = mfem.Vector(sol.Size())
        euler.Mult(sol, z)

        my_max_char_speed = euler.GetMaxCharSpeed()
        max_char_speed = MPI.COMM_WORLD.allreduce(
            my_max_char_speed, op=MPI.MAX)

        dt = cfl * hmin / max_char_speed / (2 * order + 1)

    t = 0.0
    euler.SetTime(t)
    ode_solver.Init(euler)

    # Integrate in time.
    done = False
    ti = 0
    while not done:
        dt_real = min(dt, t_final - t)
        t, dt_real = ode_solver.Step(sol, t, dt_real)

        if (cfl > 0):
            my_max_char_speed = euler.GetMaxCharSpeed()
            max_char_speed = MPI.COMM_WORLD.allreduce(
                my_max_char_speed, op=MPI.MAX)

            dt = cfl * hmin / max_char_speed / (2*order+1)

        ti = ti+1
        done = (t >= t_final - 1e-8*dt)
        if (done or ti % vis_steps == 0):
            if myid == 0:
                print("time step: " + str(ti) + ", time: " + "{:g}".format(t))
            if (visualization):
                sout << "window_title 'momentum, t = " << "{:g}".format(
                    t) << "'\n"
                sout << "parallel " << str(
                    num_procs) << " " << str(myid) << "\n"
                sout << "solution\n" << pmesh << mom
                sout.flush()

    if myid == 0:
        print("done")

    # 9. Save the final solution. This output can be viewed later using GLVis:
    #    "glvis -np 4 -m euler-mesh-final -g euler-1-final" (for x-momentum).
    pmesh.Print("euler-mesh-final."+smyid, 8)
    for k in range(num_equation):
        uk = mfem.ParGridFunction(fes, mfem.Vector(
            sol.GetDataArray()[k*fes.GetNDofs():]))
        sol_name = "euler-" + str(k) + "-final."+smyid
        uk.Save(sol_name, 8)

    # 10. Compute the L2 solution error summed for all components.
    if True:
        error = sol.ComputeLpError(2., u0)
        if myid == 0:
            print("Solution error: " + "{:g}".format(error))


if __name__ == "__main__":

    parser = ArgParser(description='Ex18p')
    parser.add_argument('-m', '--mesh',
                        default='',
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
    parser.add_argument('-novis', '--no_visualization',
                        action='store_true', default=False,
                        help='Disable GLVis visualization')
    parser.add_argument("-ea", "--element-assembly-divergence",
                        action='store_true', default=False,
                        help="Weak divergence assembly level\n" +
                        "    ea - Element assembly with interpolated")
    parser.add_argument("-mf", "--matrix-free-divergence",
                        action='store_true', default=False,
                        help="Weak divergence assembly level\n" +
                        "    mf - Nonlinear assembly in matrix-free manner")
    parser.add_argument('-vs', '--visualization-steps',
                        action='store', default=50, type=float,
                        help="Visualize every n-th timestep.")

    args = parser.parse_args()

    visualization = not args.no_visualization

    if (not args.matrix_free_divergence and
            not args.element_assembly_divergence):
        args.element_assembly_divergence = True
        args.matrix_free_divergence = False
        preassembleWeakDiv = True

    elif args.element_assembly_divergence:
        args.matrix_free_divergence = False
        preassembleWeakDiv = True

    elif args.matrix_free_divergence:
        args.element_assembly_divergence = False
        preassembleWeakDiv = False

    if myid == 0:
        parser.print_options(args)

    run(problem=args.problem,
        ser_ref_levels=args.refine_serial,
        par_ref_levels=args.refine_parallel,
        order=args.order,
        ode_solver_type=args.ode_solver,
        t_final=args.t_final,
        dt=args.time_step,
        cfl=args.cfl_number,
        visualization=visualization,
        vis_steps=args.visualization_steps,
        preassembleWeakDiv=preassembleWeakDiv,
        meshfile=args.mesh)
