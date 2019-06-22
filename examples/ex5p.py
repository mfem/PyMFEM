'''
   MFEM example 5p

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join
from mpi4py import MPI
import numpy as np
from numpy import sin, cos, exp

#### time.clock deprecated and removed in PY3.8
import time
try:
    clock = time.process_time
except AttributeError:
    clock = time.clock

def pFunc_exact(x):
    xi = float(x[0]); yi = float(x[1]); zi = 0.0
    if len(x) == 3: zi = x[2]
    from numpy import sin, cos, exp        
    return exp(xi)*sin(yi)*cos(zi)

class uFunc_ex(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       xi = float(x[0]); yi = float(x[1]); zi = 0.0
       if len(x) == 3: zi = x[2]
       ret = [- exp(xi)*sin(yi)*cos(zi),
              - exp(xi)*cos(yi)*cos(zi)]
       if len(x) == 3:
          ret.append(exp(xi)*sin(yi)*sin(zi))
       return ret
       
class pFunc_ex(mfem.PyCoefficient):
   def EvalValue(self, x):
       return  pFunc_exact(x)
   
class fFunc(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       if len(x) == 3:
           return [0., 0., 0.]
       else:
           return [0., 0.]

class gFunc(mfem.PyCoefficient):
   def EvalValue(self, x):
       if len(x) == 3:
           return  -pFunc_exact(x)
       else:
           return 0.0
class f_natural(mfem.PyCoefficient):
   def EvalValue(self, x):    
       return -pFunc_exact(x)

num_proc = MPI.COMM_WORLD.size
myid     = MPI.COMM_WORLD.rank
verbose = (myid == 0)

order =1 
static_cond = False
meshfile = expanduser(join(path, 'data', 'star.mesh'))
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()

ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
   mesh.UniformRefinement();

pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
par_ref_levels = 2
for l in range(par_ref_levels): pmesh.UniformRefinement()

hdiv_coll = mfem.RT_FECollection(order, dim)
l2_coll = mfem.L2_FECollection(order, dim)

R_space = mfem.ParFiniteElementSpace(pmesh, hdiv_coll)
W_space = mfem.ParFiniteElementSpace(pmesh, l2_coll)

dimR = R_space.GlobalTrueVSize()
dimW = W_space.GlobalTrueVSize()

if verbose:
    print("***********************************************************")
    print("dim(R) = " + str(dimR))
    print("dim(W) = " + str(dimW))
    print("dim(R+W) = " + str(dimR+dimW))
    print("***********************************************************")

block_offsets = intArray([0, R_space.GetVSize(), W_space.GetVSize()])
block_offsets.PartialSum()
block_trueOffsets = intArray([0, R_space.TrueVSize(), W_space.TrueVSize()])
block_trueOffsets.PartialSum()

k = mfem.ConstantCoefficient(1.0)

fcoeff = fFunc(dim)
fnatcoeff = f_natural()
gcoeff = gFunc()
ucoeff = uFunc_ex(dim)
pcoeff = pFunc_ex()

x = mfem.BlockVector(block_offsets)

rhs = mfem.BlockVector(block_offsets)
trueX = mfem.BlockVector(block_trueOffsets)

trueRhs = mfem.BlockVector(block_trueOffsets)
trueRhs.Assign(0.0)

fform = mfem.ParLinearForm()
fform.Update(R_space, rhs.GetBlock(0), 0)
fform.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(fcoeff))
fform.AddBoundaryIntegrator(mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
fform.Assemble();
fform.ParallelAssemble(trueRhs.GetBlock(0))

gform = mfem.ParLinearForm()
gform.Update(W_space, rhs.GetBlock(1), 0)
gform.AddDomainIntegrator(mfem.DomainLFIntegrator(gcoeff))
gform.Assemble()
gform.ParallelAssemble(trueRhs.GetBlock(1))

mVarf = mfem.ParBilinearForm(R_space)
bVarf = mfem.ParMixedBilinearForm(R_space, W_space)


mVarf.AddDomainIntegrator(mfem.VectorFEMassIntegrator(k))
mVarf.Assemble();
mVarf.Finalize();
M = mVarf.ParallelAssemble()

bVarf.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator())
bVarf.Assemble()
bVarf.Finalize()
B = bVarf.ParallelAssemble()
B *= -1;
BT = B.Transpose()

darcyOp = mfem.BlockOperator(block_trueOffsets)
darcyOp.SetBlock(0,0, M)
darcyOp.SetBlock(0,1, BT)
darcyOp.SetBlock(1,0, B)

#M2 = M.Transpose()
#M3 = M2.Transpose()
MinvBt =  B.Transpose()
Md = mfem.HypreParVector(MPI.COMM_WORLD, M.GetGlobalNumRows(),
                                         M.GetRowStarts())
M.GetDiag(Md)
MinvBt.InvScaleRows(Md)
S = mfem.hypre.ParMult(B, MinvBt)

invM = mfem.HypreDiagScale(M)
invS = mfem.HypreBoomerAMG(S)
invM.iterative_mode = False
invS.iterative_mode = False

darcyPr = mfem.BlockDiagonalPreconditioner(block_trueOffsets)
darcyPr.SetDiagonalBlock(0, invM);
darcyPr.SetDiagonalBlock(1, invS);

maxIter = 500; rtol = 1e-6; atol = 1e-10

import time
stime = clock()
solver = mfem.MINRESSolver(MPI.COMM_WORLD)

solver.SetAbsTol(atol)
solver.SetRelTol(rtol)
solver.SetMaxIter(maxIter)
solver.SetOperator(darcyOp)
solver.SetPreconditioner(darcyPr)
solver.SetPrintLevel(1)
trueX.Assign(0.0)
solver.Mult(trueRhs, trueX)

solve_time = clock() - stime

if verbose:
    if solver.GetConverged():
       print("MINRES converged in " + str(solver.GetNumIterations()) + 
          " iterations with a residual norm of " + str(solver.GetFinalNorm()))
    else:
       print("MINRES did not converge in " + str(solver.GetNumIterations()) + 
             " iterations. Residual norm is " + str(solver.GetFinalNorm()))
    print("MINRES solver took " + str(solve_time) +  "s.")

u = mfem.ParGridFunction(); p = mfem.ParGridFunction()
u.MakeRef(R_space, x.GetBlock(0), 0)
p.MakeRef(W_space, x.GetBlock(1), 0)
u.Distribute(trueX.GetBlock(0))
p.Distribute(trueX.GetBlock(1))

order_quad = max(2, 2*order+1);
irs = [mfem.IntRules.Get(i, order_quad)
       for i in range(mfem.Geometry.NumGeom)]
       
err_u  = u.ComputeL2Error(ucoeff, irs);
norm_u = mfem.ComputeGlobalLpNorm(2, ucoeff, pmesh, irs)
err_p  = p.ComputeL2Error(pcoeff, irs);
norm_p = mfem.ComputeGlobalLpNorm(2, pcoeff, pmesh, irs);

if verbose:
    print("|| u_h - u_ex || / || u_ex || = " + str(err_u / norm_u))
    print("|| p_h - p_ex || / || p_ex || = " + str(err_p / norm_p))

