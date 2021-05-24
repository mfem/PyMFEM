'''
   MFEM example 4

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname, exists
import numpy as np
from numpy import sin, array, cos

set_bc = True
static_cond = False
hybridization = False

freq = 1.0
kappa = np.pi*freq

order = 1


meshfile = expanduser(join(path, 'data', 'star.mesh'))
if not exists(meshfile):
    path = dirname(dirname(__file__))
    path = dirname(path)    
    meshfile = expanduser(join(path, 'data', 'star.mesh'))

mesh = mfem.Mesh(meshfile, 1,1)
dim = mesh.Dimension()
sdim = mesh.SpaceDimension()
ndim = mesh.Dimension()

from numba import cfunc, carray    

@cfunc(mfem.vector_sig)
def E_exact(p, out, sdim, vdim):
    out_array = carray(out, (vdim, ))    
    x = p[0]; y = p[1]
    F0 = cos(kappa*x)*sin(kappa*y)
    F1 = cos(kappa*y)*sin(kappa*x)
    out_array[0] = F0
    out_array[1] = F1
    if vdim == 3:
        out_array[2] = 0.0
        
@cfunc(mfem.vector_sig)        
def f_exact(p, out, sdim, vdim):
    out_array = carray(out, (vdim, ))    
    x = p[0]; y = p[1]
    temp = 1. + 2.*kappa*kappa
       
    F0 = temp * cos(kappa*x)*sin(kappa*y)
    F1 = temp * cos(kappa*y)*sin(kappa*x)
    out_array[0] = F0
    out_array[1] = F1
    if vdim == 3:
        out_array[2] = 0.0

ref_levels = int(np.floor(np.log(25000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
   mesh.UniformRefinement();

fec = mfem.RT_FECollection(order-1, dim)
fespace = mfem.FiniteElementSpace(mesh, fec)

print("Number of finite element unknows : " + str(fespace.GetTrueVSize()))

ess_tdof_list = intArray();
if mesh.bdr_attributes.Size():
    ess_bdr = intArray(mesh.bdr_attributes.Max())
    if set_bc: ess_bdr.Assign(1)
    else: ess_bdr.Assign(0)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

b = mfem.LinearForm(fespace);
f = mfem.VectorNumbaFunction(f_exact, sdim, ndim).GenerateCoefficient()
dd = mfem.VectorFEDomainLFIntegrator(f)
b.AddDomainIntegrator(dd)
b.Assemble()
    
x = mfem.GridFunction(fespace)
F = mfem.VectorNumbaFunction(E_exact, sdim, ndim).GenerateCoefficient()
x.ProjectCoefficient(F);

alpha = mfem.ConstantCoefficient(1.0);
beta = mfem.ConstantCoefficient(1.0);
a = mfem.BilinearForm(fespace);
a.AddDomainIntegrator(mfem.DivDivIntegrator(alpha));
a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(beta));


if (static_cond):  a.EnableStaticCondensation()
elif (hybridization):
    hfec = mfem.DG_Interface_FECollection(order-1, dim)
    hfes = mfem.FiniteElementSpace(mesh, hfec)
    a.EnableHybridization(hfes, mfem.NormalTraceJumpIntegrator(),
                          ess_tdof_list)
a.Assemble();


A = mfem.OperatorPtr()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
## Here, original version calls hegith, which is not
## defined in the header...!?
print("Size of linear system: " + str(A.Height())) 

# 10. Solve
AA = mfem.OperatorHandle2SparseMatrix(A)
M = mfem.GSSmoother(AA)
X.Print("x")
mfem.PCG(AA, M, B, X, 1, 10000, 1e-20, 0.0);

# 11. Recover the solution as a finite element grid function.
a.RecoverFEMSolution(X, b, x);

print("|| F_h - F ||_{L^2} = " + str(x.ComputeL2Error(F)))

mesh.Print('refined.mesh', 8)
x.Save('sol.gf', 8)
