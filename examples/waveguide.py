#
#   This example shows how to add waveguide port boundary 
#   condtion in plain MFEM and PyMFEM
#
#   This model simulates stragiht waveguide with 2 ports  
#   Bdry 2 : excitaiton port
#   Bdry 5 : passive port
#
#   Only TE01 modes are concerned, thought including multiple
#   higher modes is not complicated
#
#from ifigure.interactive import figure

# excitation port
ext = 'port2'

from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join
import numpy as np

epsilon0 = 8.8541878176e-12       # vacuum permittivity
mu0      = 4* np.pi*1e-7          # vacuum permiability
c        = 2.99792458e8           # speed of light
f        = 4.6e9                    # 5GHz
omega = 2*f*np.pi                 # omega =  2*pi*f
sigma = -0.0                        # resistivity
repsilon = 1.0                    # relative permittivity
rmu      = 1.0                    # relative permiability


order = 1
meshfile = expanduser(join(path, 'data', 'waveguide_hex.mesh'))
#meshfile = expanduser('~/src/PyMFEM/data/waveguide.mesh')
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()
sdim= mesh.SpaceDimension()

#   3. Refine the mesh to increase the resolution. In this example we do
#      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
#      largest number that gives a final mesh with no more than 50,000
#      elements.
#ref_levels = int(np.floor(np.log(50000./mesh.GetNE())/np.log(2.)/dim))
#for x in range(ref_levels):
for x in range(0):
     mesh.UniformRefinement();
mesh.ReorientTetMesh();


#  4. Define a finite element space on the mesh. Here we use the Nedelec
#     finite elements of the specified order.

fec = mfem.ND_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fec)

print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))

# 5. Determine the list of true (i.e. conforming) essential boundary dofs.
#    In this example, the boundary conditions are defined by marking all
#    the boundary attributes from the mesh as essential (Dirichlet) and
#    converting them to a list of true dofs.

ess_tdof_list = intArray();
if mesh.bdr_attributes.Size():
    ess_bdr = intArray([1, 0, 1, 1, 0, 1])
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
# 6. Set up the linear form b(.) which corresponds to the right-hand side
#    of the FEM linear system, which in this case is (f,phi_i) where f is
#    given by the function f_exact and phi_i are the basis functions in the
#    finite element fespace.

#from mfem.examples.waveguide import r_Jt_port_cb, i_Jt_port_cb, Ht_port_cb
r_b = mfem.LinearForm(fespace);
i_b = mfem.LinearForm(fespace);
port1_bdr = intArray([0, 1, 0, 0, 0, 0])
port2_bdr = intArray([0, 0, 0, 0, 1, 0])
#r_Jt = mfem.VectorFunctionCoefficient(sdim, r_Ht_port_cb)


m  = 0
n  = 1
a  = 0.007
b  = 0.06
k0 = 2*np.pi*f/c
k = np.sqrt(k0*k0 - np.pi*np.pi/b/b)
C = b/2
class C_Et(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi), 0, 0]
#class C_Et2(mfem.VectorPyCoefficient):
#   def EvalValue(self, x):
#       return [np.sin(2*x[2]/b*np.pi), 0, 0]
class C_Ht1(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [0., 0., k/mu0*np.sin(x[2]/b*np.pi)]
class C_Ht2(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [0., 0., -k/mu0*np.sin(x[2]/b*np.pi)]

#r_Jt = C_Jt(sdim)
if ext == 'port1':
   i_Ht = C_Ht1(sdim)
   i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port1_bdr)
else:
   i_Ht = C_Ht2(sdim)
   i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port2_bdr)
#r_dd = mfem.VectorFEBoundaryTangentLFIntegrator(r_f);
i_dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f); #\int test(E)*nxH dS
#r_b.AddBoundaryIntegrator(r_dd)
r_b.Assemble();
i_b.AddBoundaryIntegrator(i_dd)
i_b.Assemble();

active_port = intArray([0, 1, 0, 0, 0, 0])
pt1 = mfem.LinearForm(fespace);
i_Ht = C_Ht1(sdim)
i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port1_bdr)
dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f)
pt1.AddBoundaryIntegrator(dd)
pt1.Assemble()

pt2 = mfem.LinearForm(fespace);
r_Et = C_Et(sdim)
r_f =  mfem.VectorRestrictedCoefficient(r_Et, port1_bdr)
dd = mfem.VectorFEDomainLFIntegrator(r_f)
pt2.AddBoundaryIntegrator(dd)
pt2.Assemble()
r_x = mfem.GridFunction(fespace)
r_x.Assign(0.0)
r_x.ProjectBdrCoefficientTangent(r_Et, port1_bdr)
pt2d = pt2.GetDataArray()/np.sum(pt2.GetDataArray()*r_x.GetDataArray())

passive_port = intArray([0, 0, 0, 0, 1, 0])
pt3 = mfem.LinearForm(fespace);
i_Ht = C_Ht2(sdim)
i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port2_bdr)
dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f)
pt3.AddBoundaryIntegrator(dd)
pt3.Assemble()

pt4 = mfem.LinearForm(fespace);
r_Et = C_Et(sdim)
r_f =  mfem.VectorRestrictedCoefficient(r_Et, port2_bdr)
dd = mfem.VectorFEDomainLFIntegrator(r_f)
pt4.AddBoundaryIntegrator(dd)
pt4.Assemble()
r_x.Assign(0.0)
r_x.ProjectBdrCoefficientTangent(r_Et, port2_bdr)
pt4d = pt4.GetDataArray()/np.sum(pt4.GetDataArray()*r_x.GetDataArray())

# 7. Define the solution vector x as a finite element grid function
#    corresponding to fespace. Initialize x by projecting the exact
#    solution. Note that only values from the boundary edges will be used
#    when eliminating the non-homogeneous boundary condition to modify the
#    r.h.s. vector b.

r_x = mfem.GridFunction(fespace)
i_x = mfem.GridFunction(fespace)

E = mfem.VectorArrayCoefficient(3)
for i in range(3):
    E.Set(i, mfem.ConstantCoefficient(0.0))
# Don't do this since, E will try to release the same object 3 times..    
#   zero = mfem.ConstantCoefficient(0.0)    
#   E.Set(0, zero);E.Set(1, zero);E.Set(2, zero)
class C_EexactR(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi)*np.cos(x[1]*k), 0, 0]
class C_EexactI(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi)*np.sin(x[1]*k), 0, 0]
ER = C_EexactR(sdim)
EI = C_EexactI(sdim)
r_x.ProjectCoefficient(ER);
i_x.ProjectCoefficient(EI);
sol_exact = r_x.GetDataArray() + 1j*i_x.GetDataArray()
# 8. Set up the bilinear form corresponding to the EM diffusion operator
#       curl muinv curl + sigma I, by adding the curl-curl and the mass domain
#       integrators.


muinv = mfem.ConstantCoefficient(1./mu0/rmu);
mass = -(repsilon*epsilon0*omega*omega - 1j*omega*sigma)

co_rmass = mfem.ConstantCoefficient(mass.real)
co_imass = mfem.ConstantCoefficient(mass.imag)

r_a = mfem.BilinearForm(fespace);    # real part
i_a = mfem.BilinearForm(fespace);    # imag part

r_a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
r_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(co_rmass))
i_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(co_imass))

# 9. Assemble the bilinear form and the corresponding linear system,
#       applying any necessary transformations such as: eliminating boundary
#       conditions, applying conforming constraints for non-conforming AMR,
#       static condensation, etc.

static_cond = False
if (static_cond): 
      a.EnableStaticCondensation()
      b.EnableStaticCondensation()

r_a.Assemble()
i_a.Assemble()
r_A = mfem.SparseMatrix()
r_B = mfem.Vector()
r_X = mfem.Vector()
i_A = mfem.SparseMatrix()
i_B = mfem.Vector()
i_X = mfem.Vector()
r_a.FormLinearSystem(ess_tdof_list, r_x, r_b, r_A, r_X, r_B)
i_a.FormLinearSystem(ess_tdof_list, i_x, i_b, i_A, i_X, i_B)

from mfem.pi.solver.solver_utils import make_numpy_linear_system
A1, b1 = make_numpy_linear_system(r_A, r_B)
A2, b2 = make_numpy_linear_system(i_A, i_B)
print "data"
print len(A1.data), np.max(A1.data), np.min(A1.data)
A = A1+1j*A2
b = b1+1j*b2

import scipy
pt1d= pt1.GetDataArray()
pt3d= pt3.GetDataArray()
aaa = np.sum(pt2d*sol_exact)
bbb = np.sum(pt4d*sol_exact)

#v=figure();v.plot(pt1d,'r');v.plot(pt3d)
sol_size = A.shape[0]
BigA = scipy.sparse.bmat([[A, -1j*pt1d.reshape(-1,1), -1j*pt3d.reshape(-1,1)],
                          [pt2d.reshape(1,-1),  1,  0],
                          [pt4d.reshape(1,-1),  0,  1]], format='csr')
print "BigA element"
print np.max(-1j*pt1d),  np.min(-1j*pt1d)
print np.max(pt2d), np.min(pt2d)

BigExt = np.hstack((sol_exact, 0, bbb))
if ext == 'port1':
    Bigb = np.hstack((b, 1.0, 0.0))
else:
    Bigb = np.hstack((b, 0.0, 1.0))

#v = figure();v.plot((BigA * BigExt).imag)
#v = figure();v.plot((BigA * BigExt- Bigb).imag)
## Here, original version calls hegith, which is not
## defined in the header...!?
print("Size of linear system: " + str(i_A.Size())) 

# complete essential BC elimination here.
print len(BigA.data), np.max(BigA.data), np.min(BigA.data)        
print ess_tdof_list.Size()
for x in ess_tdof_list.ToList():
   BigA[x,x] = 1
   BigA[x,-2] = 0.0
   BigA[x,-1] = 0.0
   BigA[-1,x] = 0.0
   BigA[-2,x] = 0.0
   Bigb[x] = 0.0

print np.max(Bigb), np.min(Bigb)

# 10. Solve
from mfem.pi.solver.mumps_petsc.mumps_solve import mumps_solve
Bigsol = mumps_solve(BigA, Bigb)
sol = Bigsol[:sol_size]

print 'R1, R2', np.abs(BigExt[-2:]), Bigsol[-2:], np.abs(Bigsol[-2:])
r_X.SetVector(mfem.Vector(list(sol.real)), 0)
i_X.SetVector(mfem.Vector(list(sol.imag)), 0)
# 11. Recover the solution as a finite element grid function.
r_a.RecoverFEMSolution(r_X, r_b, r_x)
i_a.RecoverFEMSolution(i_X, i_b, i_x)

mesh.PrintToFile('refined.mesh', 8)
i_x.SaveToFile('isoligf', 8)
r_x.SaveToFile('rsol.gf', 8)

Exyz = [None]*3
for idx in [1,2,3]:
   ReEx = mfem.Vector()
   r_x.GetNodalValues(ReEx, idx)  #index 1, 2, 3
   ImEx = mfem.Vector()
   i_x.GetNodalValues(ImEx, idx)  #index 1, 2, 3
   Exyz[idx-1] = ReEx.GetDataArray().copy() + 1j*ImEx.GetDataArray().copy()

# to return valuds in piScope
if 'ans' in locals(): ans((fespace, Exyz, r_x, i_x))
