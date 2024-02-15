'''
   navier_solver.py

   See c++ version in the MFEM library for more detail
'''

from mfem.par import intArray, doubleArray
import mfem.par as mfem
import os
import sys
from os.path import expanduser, join
from numpy import sin, cos, exp, sqrt, zeros, abs, pi
import numpy as np
from mpi4py import MPI

num_procs = 1
myid = 0

def vel(x, t, u):

    xi = x[0]
    yi = x[1]
    
    u = np.empty(1,2)
    u[0] = pi * sin(t) * pow(sin(pi * xi), 2.0) * sin(2.0 * pi * yi)
    u[1] = -(pi * sin(t) * sin(2.0 * pi * xi) * pow(sin(pi * yi), 2.0))

    return x, u

def p(x, t):

    xi = x[0]
    yi = x[1]

    return x, cos(pi * xi) * sin(t) * sin(pi * yi)


def accel(x, t, u, kinvis):

    xi = x[0]
    yi = x[1]

    u[0] = pi * sin(t) * sin(pi * xi) * sin(pi * yi) \
          * (-1.0
             + 2.0 * pow(pi, 2.0) * sin(t) * sin(pi * xi)
             * sin(2.0 * pi * xi) * sin(pi * yi)) \
          + pi \
          * (2.0 * kinvis * pow(pi, 2.0)
             * (1.0 - 2.0 * cos(2.0 * pi * xi)) * sin(t)
             + cos(t) * pow(sin(pi * xi), 2.0)) \
          * sin(2.0 * pi * yi)

    u[1] = pi * cos(pi * yi) * sin(t) \
          * (cos(pi * xi)
             + 2.0 * kinvis * pow(pi, 2.0) * cos(pi * yi)
             * sin(2.0 * pi * xi)) \
          - pi * (cos(t) + 6.0 * kinvis * pow(pi, 2.0) * sin(t)) \
          * sin(2.0 * pi * xi) * pow(sin(pi * yi), 2.0) \
          + 4.0 * pow(pi, 3.0) * cos(pi * yi) * pow(sin(t), 2.0) \
          * pow(sin(pi * xi), 2.0) * pow(sin(pi * yi), 3.0)
    
    return x, u

   
def run(ser_ref_levels=1,
        order=5,
        kinvis=1,
        t_final = 10 * 0.25e-4,
        dt = 0.25e-4,
        pa = True,
        ni = False,
        visualization = True,
        checkres = False):
    
    mesh = mfem.Mesh("../../data/inline-quad.mesh",1,1)
    mesh.EnsureNodes()
    nodes = mesh.GetNodes()
    nodes *= 2.0
    nodes -= 1.0
    
    for i in range(ser_ref_levels):
        mesh.UniformRefinement()
        
    print("Number of Elements " + str(mesh.GetNE()))

    if MPI.ROOT:
        print("Number of elements: " + str(mesh.GetNE()))
    
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    
    navsolv = mfem.navier_solver(pmesh,order,kinvis)
    mfem.navier_solver.EnablePA(pa)
    mfem.navier_solver.EnableNI(ni)
    
    u_ic = navsolv.GetCurrentVelocity()
    u_excoeff = mfem.VectorFunctionCoefficient(pmesh.Dimension(),vel)
    
    ### need to translate this stuff

    p_excoeff = mfem.FunctionCoefficient(p)

    pmesh.bdr_attributes.Max()
    attr = 1
    navsolv.AddVelDirichletBC(vel, attr)

    pmesh.attributes.Max()
    domain_attr = 1
    navsolv.AddAccelTerm(accel, domain_attr)

    ### need to translate this stuff

    t = 0
    last_step = False

    navsolv.Setup(dt)

    err_u = 0
    err_p = 0
    u_gf = navsolv.GetCurrentVelocity()
    p_gf = navsolv.GetCurrentPressure()

    while last_step == False:
        if t + dt >= t_final - dt/2:
            last_step = True
        
        navsolv.Step(t, dt, step) #t should update in here

        u_excoeff.SetTime(t)
        p_excoeff.SetTime(t)
        err_u = u_gf.ComputeL2Error(u_excoeff)
        err_p = p_gf.ComputeL2Error(p_excoeff)

        if MPI.ROOT:
            print("%11s %11s %11s %11s\n")
            print("{:,5E} {:,5E} {:,5E} {:,5E}\n".format(t,dt,err_u,err_p))



        step = step + 1
    
    navsolv.PrintTimingData()

    if checkres:
        tol = 1E-3
        if err_u > tol or err_p > tol:
            if MPI.ROOT:
                print("Result has a larger error than expected.")

if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='navier_mms (translated from miniapps/navier/navier_mms.cpp)')

    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly in serial.")
    parser.add_argument('-o', '--order',
                        action='store', default=5, type=int,
                        help="Order (degree) of the finite elements.")
    parser.add_argument('-dt', '--time-step',
                        action='store', default=0.25e-4, type=float,
                        help="Time step.")
    parser.add_argument('-tf', '--final-time',
                        action='store', default=10 * 0.25e-4, type=float,
                        help="Final time.")
    parser.add_argument('-no-pa', '--disable-pa',
                        action='store_false',
                        help="Disable partial assembly.")
    parser.add_argument('-ni', '--enable-pa',
                        action='store_true',
                        help="Enable numerical integration rules.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-cr', '--checkresult',
                        action='store_true',
                        help="Enable or disable checking of the result. Returns -1 on failure.")
   
    
    args = parser.parse_args()
    parser.print_options(args)

    # meshfile = expanduser(
    #     join(os.path.dirname(__file__), '..', '..', 'data', args.mesh))
    # numba = (args.numba == 1)

    run(ser_ref_levels = args.refine_serial,
        order=args.order,
        kinvis=1,
        t_final=args.final_time,
        dt=args.time_step,
        pa=True,
        ni=False,
        # visualization=args.visualization,
        checkres=False)
    
