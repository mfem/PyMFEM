'''
   navier_2dfocs_correction.py

  This code is an outline of how I think it would be ideal 
  to run the mfem simulation with solver in the loop. 
  Essentially, the issue is that "navier_2dfocs_restart.py"
    reinitializes the mesh and flowsolver every time we call the solver, 
    which introduces a lot of uneccesary cost.
  Here, we instead initialize the flowsolver and mesh one time.
    We then call the "step" function to advance the simulation
    one timestep. 
    After the timestep, we write save current velocity as an npy file
    so that we can correct it with the NN. The corrected velocity is then
    assigned to the solver and we step again.
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

#set up simulation to run in parallel
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)

#input parameters
ser_ref_levels=0,
order=4,
kinvis=0.001,
t_final = 5,
dt = 1e-3,
pa = True,
ni = False,
visualization = False,
numba = True
 
#choose mesh, this needs to be lower resolution eventually
mesh = mfem.Mesh("rect-cylinder.msh")

# refine mesh in MFEM.
for i in range(ser_ref_levels):
    mesh.UniformRefinement()

if MPI.ROOT:
    print("Number of elements: " + str(mesh.GetNE()))

pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh) #parallel mesh

#Setup navier solver with partial assembly
flowsolver = mfem.navier_solver.NavierSolver(pmesh,order,kinvis)
flowsolver.EnablePA(pa)

# Set initial condition
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


# Set boundary conditions. Inlet and cylinder are Dirichlet 0,
# other boundaries are natural
attr = intArray(pmesh.bdr_attributes.Max())
attr[0] = 1 #inlet
attr[4] = 1 #cylinder
flowsolver.AddVelDirichletBC(u_excoeff, attr)

#Flowsolver setup
time = 0.0
last_step = False

flowsolver.Setup(dt)

u_gf = flowsolver.GetCurrentVelocity()
p_gf = flowsolver.GetCurrentPressure()

step = 0

# this is the name for saving velocity
# each core saves as "nav2dv.000000","nav2dv.000001", ...
sol_name_v = "nav2dv."+smyid

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

#coords gives numpy array of how mesh nodes are distributed
nodes = mfem.ParGridFunction(u_gf)
pmesh.GetNodes(nodes)
coords =  nodes.GetDataArray()
np.save('nodes'+smyid,coords)

while last_step == False:
    if time + dt >= t_final - dt/2:
        last_step = True

    #Take one timestep
    time = flowsolver.Step(time, dt, step)

    ### Apply correction ###

    #save the velocity as a pargridfunction
    u_gf.Save(sol_name_v, 8)

    #define new pargridfunction of current velocities(Not sure if necessary)
    vel_gf = mfem.ParGridFunction(pmesh,sol_name_v)

    #get numpy array of current velocities
    vels = vel_gf.GetDataArray() 

    #This is where we would add correction
    #vels = corrected_vels

    #Set corrected velocity for solver
    #u_gf.Assign(corrected_vels)
    
    ### Correction Applied ###
    
    if visualization and step % 10 == 0:
        pvdc.SetCycle(step)
        pvdc.SetTime(time)
        pvdc.Save()

    if MPI.ROOT:
            print(" "*7 + "Time" + " "*10 + "dt" )
            print(f'{time:.5e} {dt:.5e} \n')


    step = step + 1

# flowsolver.PrintTimingData()

#save final outputs
mesh_name = "nav2d-mesh."+smyid
sol_name_v = "nav2dv-final."+smyid
pmesh.Print(mesh_name, 8)
u_gf.Save(sol_name_v, 8)