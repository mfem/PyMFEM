'''
   MFEM example 20
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt


m_ = 1.0
k_ = 1.0


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
    q[0] = 0.0
    p[0] = 1.0

    # 5. Create a Mesh for visualization in phase space
    nverts = 2*(nsteps+1) if visualization else 0
    nelems = nsteps if visualization else 0

    mesh = mfem.Mesh(2, nverts, nelems, 0, 3)

    x0 = mfem.Vector(3)
    x0.Assign(0.0)
    x1 = mfem.Vector(3)
    x1.Assign(0.0)
    v = mfem.intArray(4)
    # 6. Perform time-stepping
    e_mean = 0.0

    for i in range(nsteps):
        if i == 0:
            e[0] = hamiltonian(q[0], p[0], t)
            e_mean += e[0]
            if visualization:
                x1[0] = q[0]
                x1[1] = p[0]
                x1[2] = 0.0
                mesh.AddVertex(x0)
                # These are all same.
                # mesh.AddVertex(x0.GetDataArray())
                # mesh.AddVertex(x0,GetData())
                mesh.AddVertex(x1)

        #  6b. Advance the state of the system
        t, dt = siaSolver.Step(q, p, t, dt)
        e[i+1] = hamiltonian(q[0], p[0], t)
        e_mean += e[i+1]

        #  6d. Add results to GLVis visualization
        if visualization:
            x0[2] = t
            x1[0] = q[0]
            x1[1] = p[0]
            x1[2] = t
            mesh.AddVertex(x0)
            mesh.AddVertex(x1)
            v[0] = 2*i
            v[1] = 2*(i+1)
            v[2] = 2*(i+1)+1
            v[3] = 2*i+1
            mesh.AddQuad(v)
            # this also works ;D
            # mesh.AddQuad(v.ToList())
            #mesh.AddQuad(np.array(v.ToList(), dtype=np.int32))
    #  7. Compute and display mean and standard deviation of the energy
    e_mean /= (nsteps + 1)
    e_var = 0.0
    for i in range(nsteps+1):
        e_var += (e[i] - e_mean)**2
    e_var /= (nsteps + 1)
    print("\n".join(["",
                     "Mean and standard deviation of the energy",
                     "{:g}".format(e_mean) + "\t" + "{:g}".format(sqrt(e_var))]))

    #  9. Finalize the GLVis output
    if visualization:
        mesh.FinalizeQuadMesh(1)
        fec = mfem.H1_FECollection(1, 2)
        fespace = mfem.FiniteElementSpace(mesh, fec)
        energy = mfem.GridFunction(fespace)
        energy.Assign(0.0)

        for i in range(nsteps+1):
            energy[2*i+0] = e[i]
            energy[2*i+1] = e[i]

        sock = mfem.socketstream("localhost", 19916)
        sock.precision(8)
        sock << "solution\n" << mesh << energy
        sock << "window_title 'Energy in Phase Space'\n"
        sock << "keys\n maac\n" << "axis_labels 'q' 'p' 't'\n"
        sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex20 (Sympletic ODE)')
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
