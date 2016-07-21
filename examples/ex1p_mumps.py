'''
   MFEM example 1p

   Modified to call MUMPS direct solver. Key changes are

    1) Collect RHS to head node using MPI.Gatherv.
    2) Scatter solution to all compute node using MPI.Scatter,
       before calling RecoveryFEMSolution
    3) MUMPS is called in distributed assembly matrix mode 
       (ICNTL(5) = 0, ICNTL(18) = 3)
    4) On each compute node, HyperParCSRMatrix is converted
       to IJA matrix.
'''
import mfem.par as mfem
from mpi4py import MPI
from os.path import expanduser
import numpy as np

order = 2
static_cond = False
meshfile = expanduser('~/src/mfem-3.1/data/star.mesh')
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()
num_proc = MPI.COMM_WORLD.size
myid     = MPI.COMM_WORLD.rank

ref_levels = int(np.floor(np.log(50000./mesh.GetNE())/np.log(2.)/dim))
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

import mfem.pi.solver.mumps.hypre_to_mumps as hypre_to_mumps
import mfem.pi.solver.mumps.mumps_solve as mumps_solve
nnz, irn, jcn, m = hypre_to_mumps.form_mumps_local_d_array_simple(A, 0, 0)

#if (myid == 1):
#   print 'irn', mumps_solve.i_to_list(irn, nnz)[:100]
#   print 'jcn', mumps_solve.i_to_list(jcn, nnz)[:100]
#   print 'mat', mumps_solve.d_to_list(m, nnz)[:100]
if (myid == 0):
   print("Size of linear system: " + str(x.Size()))
   print("Size of linear system: " + str(A.GetGlobalNumRows()))

s = mumps_solve.dmumps_init()
s.set_icntl(5,0)
s.set_icntl(18,3)
#s.set_icntl(28,2) this does not work..
#s.set_icntl(29,2)

# Need to collect RHS to call mumps
disps = fespace.GetGlobalTDofNumber(0)
disps = MPI.COMM_WORLD.gather(disps, root = 0)
rcounts = fespace.GetTrueVSize()
rcounts = MPI.COMM_WORLD.gather(rcounts, root = 0)
recvdata = None
senddata = [B.GetDataArray(), fespace.GetTrueVSize()]
if myid ==0: 
   recvbuf = np.empty([A.GetGlobalNumRows()], dtype=B.GetDataArray().dtype)
   recvdata = [recvbuf, rcounts, disps, MPI.DOUBLE]

MPI.COMM_WORLD.Gatherv(senddata, recvdata, root = 0)
if myid ==0: 
   s.set_n(A.GetGlobalNumRows())
#   rhs = list(recvbuf)
   s.set_rhs(mumps_solve.d_array(recvbuf))

s.set_nz_loc(nnz)
s.set_irn_loc(irn)
s.set_jcn_loc(jcn)
s.set_a_loc(m)

# No outputs
if False:
       s.set_icntl(1, -1)
       s.set_icntl(2, -1)
       s.set_icntl(3, -1)
       s.set_icntl(4,  0)

s.set_job(mumps_solve.JOB_1_2_3)
mumps_solve.dmumps_call(s)
mumps_solve.dmumps_end(s)

# Need to scatter solution
senddata = None
recvdata = np.empty([fespace.GetTrueVSize()], dtype="float64")
if myid ==0: 
   sol = np.array(mumps_solve.d_to_list(s.get_rhs(), len(rhs)), dtype="float64")
   senddata = [sol, rcounts, disps, MPI.DOUBLE]
MPI.COMM_WORLD.Scatterv(senddata, recvdata, root = 0)
X.SetVector(mfem.Vector(list(recvdata)),0)

a.RecoverFEMSolution(X, b, x)

smyid = '{:0>6d}'.format(myid)
mesh_name  =  "mesh."+smyid
sol_name   =  "sol."+smyid

pmesh.PrintToFile(mesh_name, 8)
x.SaveToFile(sol_name, 8)




