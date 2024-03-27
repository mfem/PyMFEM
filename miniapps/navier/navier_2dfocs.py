'''
   navier_2dfocs.py

   2D flow over a cylinder.

   Example run with kinematic viscosity of 0.001 and visualization:
   python navier_2dfocs.py -kinvis 0.001 -vis

   paraview files will output in file '2dfoc'
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
   
def run(ser_ref_levels=0,
        order=4,
        kinvis=0.001,
        t_final = 5,
        dt = 1e-3,
        pa = True,
        ni = False,
        visualization = False,
        numba = True):
    
    mesh = mfem.Mesh("rect-cylinder.msh")
    # mesh = mfem.Mesh("../../data/inline-quad.mesh",1,1)
    
    for i in range(ser_ref_levels):
        mesh.UniformRefinement()

    if MPI.ROOT:
        print("Number of elements: " + str(mesh.GetNE()))
    
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)

    flowsolver = mfem.navier_solver.NavierSolver(pmesh,order,kinvis)
    flowsolver.EnablePA(pa)

    u_ic = flowsolver.GetCurrentVelocity()  
    if numba:
        @mfem.jit.vector(vdim=pmesh.Dimension(),td = True, interface = 'c++')
        def u_excoeff(x,t,u):
            xi = x[0]
            yi = x[1]

            if yi <= 1e-8:
                 u[1] = 1
            else:
                 u[1] = 0

            u[0] = 0.0
            
    else:
            assert False, "numba required"

    u_ic.ProjectCoefficient(u_excoeff)

    attr = intArray(pmesh.bdr_attributes.Max())
    attr[0] = 1 #inlet
    attr[4] = 1 #cylinder
    flowsolver.AddVelDirichletBC(u_excoeff, attr)

    time = 0.0
    last_step = False

    flowsolver.Setup(dt)

    u_gf = flowsolver.GetCurrentVelocity()
    p_gf = flowsolver.GetCurrentPressure()
    
    
    step = 0

    if visualization:
         pvdc = mfem.ParaViewDataCollection("2dfoc", pmesh)
         pvdc.SetDataFormat(mfem.VTKFormat_BINARY32)
         pvdc.SetHighOrderOutput(True)
         pvdc.SetLevelsOfDetail(order)
         pvdc.SetCycle(0)
         pvdc.SetTime(time)
         pvdc.RegisterField("velocity", u_gf)
         pvdc.RegisterField("pressure", p_gf)
         pvdc.Save()

    while last_step == False:
        if time + dt >= t_final - dt/2:
            last_step = True

        
        time = flowsolver.Step(time, dt, step) #t should update in here

        if visualization and step % 10 == 0:
            pvdc.SetCycle(step)
            pvdc.SetTime(time)
            pvdc.Save()

        if MPI.ROOT:
             print(" "*7 + "Time" + " "*10 + "dt" )
             print(f'{time:.5e} {dt:.5e} \n')


        step = step + 1
    
    flowsolver.PrintTimingData()



if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='navier_mms (translated from miniapps/navier/navier_mms.cpp)')

    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=0, type=int,
                        help="Number of times to refine the mesh uniformly in serial.")
    parser.add_argument('-o', '--order',
                        action='store', default=4, type=int,
                        help="Order (degree) of the finite elements.")
    parser.add_argument('-kinvis', '--kinvis',
                        action='store', default=0.001, type=float,
                        help="Kinematic viscosity.")                    
    parser.add_argument('-dt', '--time-step',
                        action='store', default=1e-3, type=float,
                        help="Time step.")
    parser.add_argument('-tf', '--final-time',
                        action='store', default=5, type=float,
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
    parser.add_argument("-n", "--numba",
                        default=1, action='store', type=int,
                        help="Use Number compiled coefficient")
   
    args = parser.parse_args()
    parser.print_options(args)

    numba = (args.numba == 1)

    run(ser_ref_levels = args.refine_serial,
        order=args.order,
        kinvis=args.kinvis,
        t_final=args.final_time,
        dt=args.time_step,
        pa=True,
        ni=False,
        visualization=args.visualization,
        numba=numba)

