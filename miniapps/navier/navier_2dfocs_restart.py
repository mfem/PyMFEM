'''
   navier_2dfocs)restart.py

   Restarting 2D flow over a cylinder with given initial velocity field.
   (Pressure field can also be provided)

   Example run with kinematic viscosity of 0.001 and visualization:
   python navier_2dfocs_restart.py -kinvis 0.001 -vis

   paraview files will output in file '2dfoc_restart'
'''

from mfem.par import intArray, doubleArray
import mfem.par as mfem
import os
import io
import sys
from os.path import expanduser, join
from numpy import sin, cos, exp, sqrt, zeros, abs, pi
import numpy as np
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)
   
def run(order=4,
        kinvis=0.001,
        t_final = 5,
        dt = 1e-3,
        pa = True,
        ni = False,
        visualization = False,
        numba = True):
    
    # not sure how to just reload the mesh
    mesh = mfem.Mesh("rect-cylinder.msh")
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)

    #reinitialize flow solver, ideally we wouldn't do this
    flowsolver = mfem.navier_solver.NavierSolver(pmesh,order,kinvis)
    flowsolver.EnablePA(pa)

    # I think the below 2 lines are unclear in python. 
    # u_ic points to the current flowsolver velocity, which means
    # changing u_ic sets the velocity of the flowsolver.
    # E.g. u_ic.Assign(ParGridFunction(...)) will set the initial condition

    u_ic = flowsolver.GetCurrentVelocity()
    p_ic = flowsolver.GetCurrentPressure()

    #still need to define this function just to set boundary conditions later
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

    
    #Reload velocity and pressure, then assign them
    ICv = mfem.ParGridFunction(pmesh,'nav2dv.'+smyid)
    ICp = mfem.ParGridFunction(pmesh,'nav2dp.'+smyid)
    u_ic.Assign(ICv)
    p_ic.Assign(ICp)

    #Set boundary conditions again
    attr = intArray(pmesh.bdr_attributes.Max())
    attr[0] = 1 #inlet
    attr[4] = 1 #cylinder
    flowsolver.AddVelDirichletBC(u_excoeff, attr)

    #everything below here is the same
    time = 0
    last_step = False

    flowsolver.Setup(dt)

    u_gf = flowsolver.GetCurrentVelocity()
    p_gf = flowsolver.GetCurrentPressure()
    
    step = 0

    if visualization:
         pvdc = mfem.ParaViewDataCollection("2dfoc_restart", pmesh)
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
    
    # flowsolver.PrintTimingData()

    # save mesh, velocity, and pressure
    mesh_name = "nav2d-mesh."+smyid
    sol_name_v = "nav2dv."+smyid
    sol_name_p = "nav2dp."+smyid
    pmesh.Print(mesh_name, 8)
    u_gf.Save(sol_name_v, 8)
    p_gf.Save(sol_name_p, 8)




if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='navier_2dfocs')

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

    run(order=args.order,
        kinvis=args.kinvis,
        t_final=args.final_time,
        dt=args.time_step,
        pa=True,
        ni=False,
        visualization=args.visualization,
        numba=numba)

