'''
   MFEM example 13p

    This example code solves the Maxwell (electromagnetic)
    eigenvalue problem curl curl E = lambda E with homogeneous
    Dirichlet boundary conditions E x n = 0.

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mfem.par import intArray, doubleArray
from os.path import expanduser, join
from mpi4py import MPI
import numpy as np
from numpy import sin, cos, exp

ser_ref_levels = 2
par_ref_levels = 1
order = 1
nev = 5
visualization = 1

num_proc = MPI.COMM_WORLD.size
myid     = MPI.COMM_WORLD.rank

meshfile = expanduser(join(path, 'data', 'beam-tet.mesh'))
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()

for k in range(ser_ref_levels):
    mesh.UniformRefinement()

pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
for k in range(par_ref_levels):
    pmesh.UniformRefinement()
pmesh.ReorientTetMesh()
fec = mfem.ND_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)
size = fespace.GlobalTrueVSize()
if (myid == 0):
    print("Number of unknowns: " + str(size))

one = mfem.ConstantCoefficient(1.0);
ess_bdr = intArray()
if (pmesh.bdr_attributes.Size() > 0):
   ess_bdr.SetSize(pmesh.bdr_attributes.Max())
   ess_bdr.Assign(1)

a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.CurlCurlIntegrator(one))

if (pmesh.bdr_attributes.Size() == 0):
     # Add a mass term if the mesh has no boundary, e.g. periodic mesh or
     # closed surface.
     a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(one))

a.Assemble();
a.EliminateEssentialBCDiag(ess_bdr, 1.0);
a.Finalize();

m = mfem.ParBilinearForm(fespace)
m.AddDomainIntegrator(mfem.VectorFEMassIntegrator(one))
m.Assemble()
# shift the eigenvalue corresponding to eliminated dofs to a large value
m.EliminateEssentialBCDiag(ess_bdr, 2.3e-308)
m.Finalize();

A = a.ParallelAssemble()
M = m.ParallelAssemble()
   
ams = mfem.HypreAMS(A,fespace)
ams.SetPrintLevel(0);
ams.SetSingularProblem();

ame = mfem.HypreAME(MPI.COMM_WORLD)
ame.SetNumModes(nev)
ame.SetPreconditioner(ams)
ame.SetMaxIter(100)
ame.SetTol(1e-8)
ame.SetPrintLevel(1)
ame.SetMassMatrix(M)
ame.SetOperator(A)

eigenvalues = doubleArray()
ame.Solve();
ame.GetEigenvalues(eigenvalues);
x = mfem.ParGridFunction(fespace)


smyid = '{:0>6d}'.format(myid)
mesh_name  =  "ex13_mesh."+smyid


pmesh.Print(mesh_name, 8)

for i in range(nev):
    if ( myid == 0 ):
        print("Eigenmode " + str(i+1)  +'/' + str(nev) +
              ", Lambda = " + str(eigenvalues[i]))
    
    sol_name   =  "ex13_mode_"+str(i)+"."+smyid    
    x.Assign(ame.GetEigenvector(i))
    x.Save(sol_name, 8)    
    c  = None
    if (myid == 0):
        from builtins import input
        c = input("press (q)uit or (c)ontinue --> ")
    c = MPI.COMM_WORLD.bcast(c, root=0)
    if (c != 'c'): break


