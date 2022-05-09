'''
   MFEM example 14

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np

ref_levels = -1
order = 1
sigma = -1.0
kappa = -1.0

if (kappa < 0):
    kappa = (order+1)**2.

meshfile = expanduser(
    join(dirname(__file__), '..', 'data', 'star.mesh'))
mesh = mfem.Mesh(meshfile, 1, 1)

dim = mesh.Dimension()

#   3. Refine the mesh to increase the resolution. In this example we do
#      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
#      largest number that gives a final mesh with no more than 50,000
#      elements.
if ref_levels < 0:
    ref_levels = int(np.floor(np.log(50000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

if (mesh.NURBSext):
    mesh.SetCurvature(max(order, 1))

# 4. Define a finite element space on the mesh. Here we use discontinuous
#    finite elements of the specified order >= 0.
fec = mfem.DG_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fec)
print('Number of finite element unknowns: ' + str(fespace.GetVSize()))

# 5. Set up the linear form b(.) which corresponds to the right-hand side of
#    the FEM linear system.
b = mfem.LinearForm(fespace)
one = mfem.ConstantCoefficient(1.0)
zero = mfem.ConstantCoefficient(0.0)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
b.AddBdrFaceIntegrator(
    mfem.DGDirichletLFIntegrator(zero, one, sigma, kappa))
b.Assemble()

# 6. Define the solution vector x as a finite element grid function
#    corresponding to fespace. Initialize x with initial guess of zero.
x = mfem.GridFunction(fespace)
x.Assign(0.0)

# 7. Set up the bilinear form a(.,.) on the finite element space
#    corresponding to the Laplacian operator -Delta, by adding the Diffusion
#    domain integrator and the interior and boundary DG face integrators.
#    Note that boundary conditions are imposed weakly in the form, so there
#    is no need for dof elimination. After assembly and finalizing we
#    extract the corresponding sparse matrix A.
a = mfem.BilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))

a.AddInteriorFaceIntegrator(mfem.DGDiffusionIntegrator(one, sigma, kappa))
a.AddBdrFaceIntegrator(mfem.DGDiffusionIntegrator(one, sigma, kappa))
a.Assemble()
a.Finalize()
A = a.SpMat()

# 8. Define a simple symmetric Gauss-Seidel preconditioner and use it to
#    solve the system Ax=b with PCG in the symmetric case, and GMRES in the
#    non-symmetric one.
M = mfem.GSSmoother(A)
if (sigma == -1.0):
    mfem.PCG(A, M, b, x, 1, 500, 1e-12, 0.0)
else:
    mfem.GMRES(A, M, b, x, 1, 500, 10, 1e-12, 0.0)

# 9. Save the refined mesh and the solution. This output can be viewed later
#    using GLVis: "glvis -m refined.mesh -g sol.gf".
mesh.Print('refined.mesh', 8)
x.Save('sol.gf', 8)
