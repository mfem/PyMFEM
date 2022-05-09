'''
   MFEM example 20p
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)

m_ = 1.0
k_ = 1.0


def nicePrint(*s):
    MPI.COMM_WORLD.Barrier()
    for i in range(num_procs):
        MPI.COMM_WORLD.Barrier()
        if i == myid:
            print(str(myid)+': ' + ': '.join([str(ss) for ss in s]))
        MPI.COMM_WORLD.Barrier()


def run(order=1,
        prob=0,
        nsteps=100,
        dt=0.1,
        sc=1.0,
        visualization=False):

    class GradT(mfem.Operator):
        def __init__(self):
            mfem.Operator.__init__(self, 1)

        def Mult(self, x, y):
            y.Set(1.0/m_, x)

    class NegGradV(mfem.TimeDependentOperator):
        def __init__(self):
            mfem.TimeDependentOperator.__init__(self, 1)

        def Mult(self, x, y):
            if prob == 1:
                y[0] = - k_ * sin(x[0])
            elif prob == 2:
                y[0] = - k_ * x[0] * exp(-0.5 * x[0] * x[0])
            elif prob == 3:
                y[0] = - k_ * (1.0 + 2.0 * x[0] * x[0]) * x[0]
            elif prob == 4:
                y[0] = - k_ * (1.0 - 0.25 * x[0] * x[0]) * x[0]
            else:
                y[0] = - k_ * x[0]

    def hamiltonian(q, p, t):
        h = 1.0 - 0.5 / m_ + 0.5 * p * p / m_
        if prob == 1:
            h += k_ * (1.0 - cos(q))
        elif prob == 2:
            h += k_ * (1.0 - exp(-0.5 * q * q))
        elif prob == 3:
            h += 0.5 * k_ * (1.0 + q * q) * q * q
        elif prob == 4:
            h += 0.5 * k_ * (1.0 - 0.125 * q * q) * q * q
        else:
            h += 0.5 * k_ * q * q
        return h
    
    device = mfem.Device('cpu')
    if myid == 0:        
        device.Print()

    # 2. Create and Initialize the Symplectic Integration Solver
    siaSolver = mfem.SIAVSolver(order)
    P = GradT()
    F = NegGradV()
    siaSolver.Init(P, F)

    # 3. Set the initial conditions
    t = 0.0
    q = mfem.Vector(1)
    p = mfem.Vector(1)
    e = mfem.Vector(nsteps+1)
    q[0] = sin(2*pi*myid/num_procs)
    p[0] = cos(2*pi*myid/num_procs)

    # 5. Create a Mesh for visualization in phase space
    nverts = 2*(nsteps+1)*num_procs if visualization else 0
    nelems = nsteps*num_procs if visualization else 0

    mesh = mfem.Mesh(2, nverts, nelems, 0, 3)

    part = mfem.intArray(nelems)

    # 6. Perform time-stepping
    e_mean = 0.0

    for i in range(nsteps):
        if i == 0:
            e[0] = hamiltonian(q[0], p[0], t)
            e_mean += e[0]
            if visualization:
                for j in range(num_procs):
                    mesh.AddVertex([0, 0, 0])
                    mesh.AddVertex([q[0], p[0], 0.0])

        #  6b. Advance the state of the system
        t, dt = siaSolver.Step(q, p, t, dt)
        e[i+1] = hamiltonian(q[0], p[0], t)
        e_mean += e[i+1]

        #  6d. Add results to GLVis visualization
        if visualization:
            for j in range(num_procs):
                mesh.AddVertex([0,  0, t])
                mesh.AddVertex([q[0], p[0], t])
                mesh.AddQuad([2*i*num_procs + 2*j,
                              2*(i+1)*num_procs + 2*j,
                              2*(i+1)*num_procs + 2*j+1,
                              2*i*num_procs + 2*j+1])
                part[num_procs*i + j] = j
            # this also works ;D
            # mesh.AddQuad(v.ToList())
            #mesh.AddQuad(np.array(v.ToList(), dtype=np.int32))
    #  7. Compute and display mean and standard deviation of the energy
    e_mean /= (nsteps + 1)
    e_var = 0.0
    for i in range(nsteps+1):
        e_var += (e[i] - e_mean)**2
    e_var /= (nsteps + 1)
    if myid == 0:
        print("Mean and standard deviation of the energy")
    nicePrint("{:g}".format(e_mean) + "\t" + "{:g}".format(sqrt(e_var)))

    #  9. Finalize the GLVis output
    if visualization:
        mesh.FinalizeQuadMesh(1)

        pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh, part.GetData())

        fec = mfem.H1_FECollection(1, 2)
        fespace = mfem.ParFiniteElementSpace(pmesh, fec)
        energy = mfem.ParGridFunction(fespace)
        energy.Assign(0.0)

        for i in range(nsteps+1):
            energy[2*i+0] = e[i]
            energy[2*i+1] = e[i]

        sock = mfem.socketstream("localhost", 19916)
        sock.precision(8)
        sock << "parallel " << num_procs << " " << myid << "\n"
        sock << "solution\n" << pmesh << energy
        sock << "window_title 'Energy in Phase Space'\n"
        sock << "keys\n maac\n" << "axis_labels 'q' 'p' 't'\n"
        sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex20p (Sympletic ODE)')
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-p",
                        "--problem-type",
                        action='store', type=int, default=0,
                        help=''.join(["Problem Type:\n",
                                      "\t  0 - Simple Harmonic Oscillator\n",
                                      "\t  1 - Pendulum\n",
                                      "\t  2 - Gaussian Potential Well\n",
                                      "\t  3 - Quartic Potential\n",
                                      "\t  4 - Negative Quartic Potential", ]))
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Time integration order")
    parser.add_argument('-n', '--number-of-steps',
                        action='store', default=100, type=int,
                        help="Number of time steps")
    parser.add_argument('-dt', '--time-step',
                        action='store', default=0.1, type=float,
                        help="Time step size")
    parser.add_argument('-k', '--spring-constant',
                        action='store', default=1, type=float,
                        help="Sprint constant")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=True,
                        help='Enable GLVis visualization')
    parser.add_argument('-no-gp', '--no-gnuplot',
                        action='store_true',
                        default=True,
                        help='Disable GnuPlot visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    prob = args.problem_type
    visualization = args.visualization
    order = args.order
    nsteps = args.number_of_steps
    dt = args.time_step
    sc = args.spring_constant
    np_gp = args.no_gnuplot

    run(order=order,
        prob=prob,
        nsteps=nsteps,
        dt=dt,
        sc=sc,
        visualization=visualization)
