'''
   MFEM example 1p

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mpi4py import MPI
from os.path import expanduser, join
import numpy as np

order = 1
static_cond = False
meshfile = expanduser(join(path, 'data', 'star.mesh'))
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()
num_proc = MPI.COMM_WORLD.size
myid     = MPI.COMM_WORLD.rank

ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
   mesh.UniformRefinement();
mesh.ReorientTetMesh();
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh

par_ref_levels = 2
for l in range(par_ref_levels):
    pmesh.UniformRefinement();

if order > 0:
    fec = mfem.H1_FECollection(order, dim)
elif mesh.GetNodes():
    fec = mesh.GetNodes().OwnFEC()
    prinr( "Using isoparametric FEs: " + str(fec.Name()));
else:
    order = 1
    fec = mfem.H1_FECollection(order, dim)

fespace =mfem.ParFiniteElementSpace(pmesh, fec)
fe_size = fespace.GlobalTrueVSize()

if (myid == 0):
   print('Number of finite element unknowns: '+  str(fe_size))

ess_tdof_list = mfem.intArray()
if pmesh.bdr_attributes.Size()>0:
    ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

#   the basis functions in the finite element fespace.
b = mfem.ParLinearForm(fespace)
one = mfem.ConstantCoefficient(1.0)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
b.Assemble();

x = mfem.ParGridFunction(fespace);
x.Assign(0.0)

a = mfem.ParBilinearForm(fespace);
a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))

if static_cond: a.EnableStaticCondensation()
a.Assemble();

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

if (myid == 0):
   print("Size of linear system: " + str(x.Size()))
   print("Size of linear system: " + str(A.GetGlobalNumRows()))

amg = mfem.HypreBoomerAMG(A)
pcg = mfem.HyprePCG(A)
pcg.SetTol(1e-12)
pcg.SetMaxIter(200)
pcg.SetPrintLevel(2)
pcg.SetPreconditioner(amg)
pcg.Mult(B, X);


a.RecoverFEMSolution(X, b, x)

smyid = '{:0>6d}'.format(myid)
mesh_name  =  "mesh."+smyid
sol_name   =  "sol."+smyid

pmesh.PrintToFile(mesh_name, 8)
x.SaveToFile(sol_name, 8)




