#
#   This example shows TEM waveguide port and
#   periodic boundary in plain MFEM and PyMFEM
#
#   This model simulates stragiht waveguide with 2 ports  
#   Periodic 3-4 8-9, and 15-16
#   Bdry 6/7 : TEM 
#   Internal boundary 10 and 17: natural
#
#   Domain 1: big Box
#   
#
from ifigure.interactive import figure
import mfem
from mfem import intArray
from os.path import expanduser
import numpy as np

epsilon0 = 8.8541878176e-12       # vacuum permittivity
mu0      = 4* np.pi*1e-7          # vacuum permiability
c        = 2.99792458e8           # speed of light
f        = 4.6e9                    # 5GHz
omega = 2*f*np.pi                 # omega =  2*pi*f
sigma = 0.                        # resistivity
epsilonr = 1.0                    # relative permittivity
epsilonr = 1.0                    # relative permittivity
mur      = 1.0                    # relative permiability


iport1 = 6
iport2 = 13

order = 1
meshfile = expanduser('~/PyMFEM/data/grill.mesh')
mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()
sdim= mesh.SpaceDimension()

#   3. Refine the mesh to increase the resolution. In this example we do
#      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
#      largest number that gives a final mesh with no more than 50,000
#      elements.
#ref_levels = int(np.floor(np.log(50000./mesh.GetNE())/np.log(2.)/dim))
#for x in range(ref_levels):
#for x in range(1):
#     mesh.UniformRefinement();
mesh.ReorientTetMesh();

#  4. Define a finite element space on the mesh. Here we use the Nedelec
#     finite elements of the specified order.

fec = mfem.ND_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fec)
print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))
#ans((fespace, mesh))
# 5. Determine the list of true (i.e. conforming) essential boundary dofs

ess_tdof_list = intArray();
ess_bdr = intArray([1]*(mesh.bdr_attributes.Max()))
for x in [3, 4, 8, 9, 15, 16, 10 ,17, 1, 20]: ess_bdr[x-1] = 0
#for x in [3, 4, 8, 9, 15, 16, ]: ess_bdr[x-1] = 0
ess_bdr[iport1] = 0
ess_bdr[iport2] = 0

if mesh.bdr_attributes.Size():
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

# 6. RHS and arrays for ports
#    Convetion is exp(ikx - iwt)
#    We set Ex = 1, 1/mu curl E = ik/mu = -i omega H
class C_Et(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [1, 0, 0]
class C_curlE_mu(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [0., 0., omega*np.sqrt(epsilon0*epsilonr/mu0/mur)]
# RHS
r_b = mfem.LinearForm(fespace)
i_b = mfem.LinearForm(fespace)

port1 = intArray([0]*(mesh.bdr_attributes.Max()))
port2 = intArray([0]*(mesh.bdr_attributes.Max()))
port1[iport1] = 1
port2[iport2] = 1
i_Ht = C_curlE_mu(sdim)
i_f1 =  mfem.VectorRestrictedCoefficient(i_Ht, port1)
i_dd1 = mfem.VectorFEBoundaryTangentLFIntegrator(i_f1)
i_f2 =  mfem.VectorRestrictedCoefficient(i_Ht, port2)
i_dd2 = mfem.VectorFEBoundaryTangentLFIntegrator(i_f2)

#r_b.AddBoundaryIntegrator(i_dd2) # for 180 deg phasing...
r_b.Assemble();
i_b.AddBoundaryIntegrator(i_dd1)
i_b.AddBoundaryIntegrator(i_dd2)
i_b.Assemble();

pt1 = mfem.LinearForm(fespace)
i_Ht = C_curlE_mu(sdim)
i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port1)
dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f)
pt1.AddBoundaryIntegrator(dd)
pt1.Assemble()

pt2 = mfem.LinearForm(fespace);
r_Et = C_Et(sdim)
r_f =  mfem.VectorRestrictedCoefficient(r_Et, port1)
dd = mfem.VectorFEDomainLFIntegrator(r_f)
pt2.AddBoundaryIntegrator(dd)
pt2.Assemble()
r_x = mfem.GridFunction(fespace)
r_x.Assign(0.0)
r_x.ProjectBdrCoefficientTangent(r_Et, port1)
pt2d = pt2.GetDataArray()/np.sum(pt2.GetDataArray()*r_x.GetDataArray())

pt3 = mfem.LinearForm(fespace);
i_Ht = C_curlE_mu(sdim)
i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port2)
dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f)
pt3.AddBoundaryIntegrator(dd)
pt3.Assemble()

pt4 = mfem.LinearForm(fespace);
r_Et = C_Et(sdim)
r_f =  mfem.VectorRestrictedCoefficient(r_Et, port2)
dd = mfem.VectorFEDomainLFIntegrator(r_f)
pt4.AddBoundaryIntegrator(dd)
pt4.Assemble()
r_x.Assign(0.0)
r_x.ProjectBdrCoefficientTangent(r_Et, port2)
pt4d = pt4.GetDataArray()/np.sum(pt4.GetDataArray()*r_x.GetDataArray())

# 7. Define the solution vector x as a finite element grid function
#    corresponding to fespace. Initialize x by projecting the exact
#    solution. Note that only values from the boundary edges will be used
#    when eliminating the non-homogeneous boundary condition to modify the
#    r.h.s. vector b.

r_x = mfem.GridFunction(fespace); r_x.Assign(0.0)
i_x = mfem.GridFunction(fespace); i_x.Assign(0.0)
sol_exact = r_x.GetDataArray() + 1j*i_x.GetDataArray()
#zero = mfem.ConstantCoefficient(0.0)
#E = mfem.VectorArrayCoefficient(3)
#E.Set(0, zero);E.Set(1, zero);E.Set(2, zero)
#r_x.ProjectCoefficient(E);
#i_x.ProjectCoefficient(E);
#sol_exact = r_x.GetDataArray() + 1j*i_x.GetDataArray()
# 8. Set up the bilinear form corresponding to the EM diffusion operator
#       curl muinv curl + sigma I, by adding the curl-curl and the mass domain
#       integrators.
#ans((fespace, mesh, r_x))
muinv = mfem.ConstantCoefficient(1./mu0/mur);

mass = -(epsilonr*epsilon0*omega*omega + 1j*omega*sigma)
mass2 = -1j*omega*0.03

co_rmass = mfem.ConstantCoefficient(mass.real)
co_imass = mfem.ConstantCoefficient(mass.imag)
ab_imass0 = mfem.ConstantCoefficient(mass2.imag)


r_a = mfem.BilinearForm(fespace);    # real part
i_a = mfem.BilinearForm(fespace);    # imag part

domains = intArray([0]*(mesh.attributes.Max()))
domains[0] = 1
ab_imass =  mfem.RestrictedCoefficient(ab_imass0, domains)

r_a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
r_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(co_rmass))
i_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(co_imass))
i_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(ab_imass))

#
# PBC (periodic boundary condition) mapping
# put -1 to diagnal elements of non-empty row 
# then eliminate empty row
from scipy.sparse import csr_matrix
map = proj.find_dof_maps(fespace, [3, 8, 15], [4, 9, 16], 
      trans = 'xy', eps= 1e-10)
map = csr_matrix(map)
num_nonzeros = np.diff(map.indptr)
for i in  np.where(num_nonzeros != 0)[0]: map[i, i] = -1
map = map[num_nonzeros != 0]

map2 = proj.find_dof_maps(fespace, [1], [20], 
      trans = 'yz', eps= 1e-10)
map2 = csr_matrix(map2, dtype=complex)
num_nonzeros = np.diff(map2.indptr)
for i in  np.where(num_nonzeros != 0)[0]: map2[i, i] = -np.exp(-1j*np.pi/180*40.)
map2 = map2[num_nonzeros != 0]

#ans((fespace, map))
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

from mfem.helpers.solver.mumps_solve import mumps_solve, make_numpy_linear_system
A1, b1 = make_numpy_linear_system(r_A, r_b)
A2, b2 = make_numpy_linear_system(i_A, i_b)

A = A1+1j*A2
b = b1+1j*b2

import scipy
pt1d= pt1.GetDataArray()
pt3d= pt3.GetDataArray()

#v=figure();v.plot(pt1d,'r');v.plot(pt3d)
sol_size = A.shape[0]
BigA = scipy.sparse.bmat([[A, -1j*pt1d.reshape(-1,1), -1j*pt3d.reshape(-1,1)],
                          [pt2d.reshape(1,-1),  1,  0],
                          [pt4d.reshape(1,-1),  0,  1]], format='csr')

#BigExt = np.hstack((sol_exact, 0, bbb))
Bigb = np.hstack((b, 1, 1))

#v = figure();v.plot((BigA * BigExt).imag)
#v = figure();v.plot((BigA * BigExt- Bigb).imag)
## Here, original version calls hegith, which is not
## defined in the header...!?
print("Size of linear system: " + str(i_A.Size())) 
print("Size of linear system: " + str(BigA.shape)) 

# 10. Solve
#ans((BigA, Bigb, mumps_solve(BigA, Bigb, c = map, Acbc = True)))
Bigsol = mumps_solve(BigA, Bigb, c = [map, map2], use_null = False)

sol = Bigsol[:sol_size]
print 'R1, R2', np.abs(Bigsol[sol_size:sol_size+2])
r_X.SetVector(mfem.Vector(list(sol.real)), 0)
i_X.SetVector(mfem.Vector(list(sol.imag)), 0)
# 11. Recover the solution as a finite element grid function.
r_a.RecoverFEMSolution(r_X, r_b, r_x)
i_a.RecoverFEMSolution(i_X, i_b, i_x)

mesh.PrintToFile('refined.mesh', 8)
i_x.SaveToFile('soli.gf', 8)
r_x.SaveToFile('solr.gf', 8)

Exyz = [None]*3
for idx in [1,2,3]:
   ReEx = mfem.Vector()
   r_x.GetNodalValues(ReEx, idx)  #index 1, 2, 3
   ImEx = mfem.Vector()
   i_x.GetNodalValues(ImEx, idx)  #index 1, 2, 3
   Exyz[idx-1] = ReEx.GetDataArray().copy() + 1j*ImEx.GetDataArray().copy()
ans((fespace, Exyz, r_x, i_x))
