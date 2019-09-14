'''
   MFEM example 5

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join
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
       if len(x) == 3: return  -pFunc_exact(x)       
       return 0.
   
class f_natural(mfem.PyCoefficient):
   def EvalValue(self, x):
       return -pFunc_exact(x)


order =1 
static_cond = False
meshfile = expanduser(join(path, 'data', 'star.mesh'))
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()

ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
   mesh.UniformRefinement();

hdiv_coll = mfem.RT_FECollection(order, dim)
l2_coll = mfem.L2_FECollection(order, dim)

R_space = mfem.FiniteElementSpace(mesh, hdiv_coll)
W_space = mfem.FiniteElementSpace(mesh, l2_coll)

dimR = R_space.GetVSize()
dimW = W_space.GetVSize()

print("***********************************************************")
print("dim(R) = " + str(dimR))
print("dim(W) = " + str(dimW))
print("dim(R+W) = " + str(dimR+dimW))
print("***********************************************************")

block_offsets = intArray([0, dimR, dimW])
block_offsets.PartialSum()

k = mfem.ConstantCoefficient(1.0)

fcoeff = fFunc(dim)
fnatcoeff = f_natural()
gcoeff = gFunc()
ucoeff = uFunc_ex(dim)
pcoeff = pFunc_ex()

x = mfem.BlockVector(block_offsets)
rhs = mfem.BlockVector(block_offsets)

fform = mfem.LinearForm()
fform.Update(R_space, rhs.GetBlock(0), 0)
fform.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(fcoeff))
fform.AddBoundaryIntegrator(mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
fform.Assemble();

gform = mfem.LinearForm()
gform.Update(W_space, rhs.GetBlock(1), 0)
gform.AddDomainIntegrator(mfem.DomainLFIntegrator(gcoeff))
gform.Assemble()

mVarf = mfem.BilinearForm(R_space)
bVarf = mfem.MixedBilinearForm(R_space, W_space)

mVarf.AddDomainIntegrator(mfem.VectorFEMassIntegrator(k))
mVarf.Assemble();
mVarf.Finalize();
M = mVarf.SpMat()

bVarf.AddDomainIntegrator( mfem.VectorFEDivergenceIntegrator())
bVarf.Assemble()
bVarf.Finalize()
B = bVarf.SpMat()
B *= -1;
BT = mfem.Transpose(B)

darcyOp = mfem.BlockOperator(block_offsets)
darcyOp.SetBlock(0,0, M)
darcyOp.SetBlock(0,1, BT)
darcyOp.SetBlock(1,0, B)

MinvBt = mfem.Transpose(B)
Md = mfem.Vector(M.Height())
M.GetDiag(Md)
for i in range(Md.Size()):
   MinvBt.ScaleRow(i, 1/Md[i])
S = mfem.Mult(B, MinvBt)

invM = mfem.DSmoother(M)
invS = mfem.GSSmoother(S)
invM.iterative_mode = False;
invS.iterative_mode = False;

darcyPrec = mfem.BlockDiagonalPreconditioner(block_offsets)
darcyPrec.SetDiagonalBlock(0, invM);
darcyPrec.SetDiagonalBlock(1, invS);

maxIter = 500; rtol = 1e-6; atol = 1e-10

stime = clock()
solver = mfem.MINRESSolver()
solver.SetAbsTol(atol)
solver.SetRelTol(rtol)
solver.SetMaxIter(maxIter)
solver.SetOperator(darcyOp)
solver.SetPreconditioner(darcyPrec)
solver.SetPrintLevel(1)
x.Assign(0.0)
solver.Mult(rhs, x)

solve_time = clock() - stime

if solver.GetConverged():
   print("MINRES converged in " + str(solver.GetNumIterations()) + 
      " iterations with a residual norm of " + str(solver.GetFinalNorm()))
else:
   print("MINRES did not converge in " + str(solver.GetNumIterations()) + 
         " iterations. Residual norm is " + str(solver.GetFinalNorm()))
print("MINRES solver took " + str(solve_time) +  "s.")

u = mfem.GridFunction(); p = mfem.GridFunction()
u.MakeRef(R_space, x.GetBlock(0), 0)
p.MakeRef(W_space, x.GetBlock(1), 0)

order_quad = max(2, 2*order+1);

irs =[mfem.IntRules.Get(i, order_quad)
      for i in range(mfem.Geometry.NumGeom)]

norm_p = mfem.ComputeLpNorm(2, pcoeff, mesh, irs)
norm_u = mfem.ComputeLpNorm(2, ucoeff, mesh, irs)
err_u  = u.ComputeL2Error(ucoeff, irs)
err_p  = p.ComputeL2Error(pcoeff, irs)

print("|| u_h - u_ex || / || u_ex || = " + str(err_u / norm_u))
print("|| p_h - p_ex || / || p_ex || = " + str(err_p / norm_p))




