'''
   MFEM example 3

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join
import numpy as np
from numpy import sin, array

freq = 1.0
kappa = np.pi*freq

static_cond = False
order = 1
meshfile = expanduser(join(path, 'data', 'beam-tet.mesh'))

mesh = mfem.Mesh(meshfile, 1,1)

dim = mesh.Dimension()
sdim= mesh.SpaceDimension()

class E_exact(mfem.VectorPyCoefficient):
    def __init__(self):
       
       mfem.VectorPyCoefficient.__init__(self, dim)
    def EvalValue(self, x):
       return (sin(kappa * x[1]),
               sin(kappa * x[2]),
               sin(kappa * x[0]))
class f_exact(mfem.VectorPyCoefficient):
    def __init__(self):
       mfem.VectorPyCoefficient.__init__(self, dim)
    def EvalValue(self, x):   
       return ((1 + kappa**2)*sin(kappa * x[1]),
               (1 + kappa**2)*sin(kappa * x[2]),
               (1 + kappa**2)*sin(kappa * x[0]))

#   3. Refine the mesh to increase the resolution. In this example we do
#      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
#      largest number that gives a final mesh with no more than 50,000
#      elements.

ref_levels = int(np.floor(np.log(50000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
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
    ess_bdr = intArray(mesh.bdr_attributes.Max())
    ess_bdr = intArray([1]*mesh.bdr_attributes.Max())
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

# 6. Set up the linear form b(.) which corresponds to the right-hand side
#    of the FEM linear system, which in this case is (f,phi_i) where f is
#    given by the function f_exact and phi_i are the basis functions in the
#    finite element fespace.

b = mfem.LinearForm(fespace);
f = f_exact()
dd = mfem.VectorFEDomainLFIntegrator(f);
b.AddDomainIntegrator(dd)
b.Assemble();

# 7. Define the solution vector x as a finite element grid function
#    corresponding to fespace. Initialize x by projecting the exact
#    solution. Note that only values from the boundary edges will be used
#    when eliminating the non-homogeneous boundary condition to modify the
#    r.h.s. vector b.

#from mfem.examples.ex3 import E_exact_cb
x = mfem.GridFunction(fespace)
E = E_exact()
x.ProjectCoefficient(E);

# 8. Set up the bilinear form corresponding to the EM diffusion operator
#       curl muinv curl + sigma I, by adding the curl-curl and the mass domain
#       integrators.

muinv = mfem.ConstantCoefficient(1.0);
sigma = mfem.ConstantCoefficient(1.0);
a = mfem.BilinearForm(fespace);
a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv));
a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma));

# 9. Assemble the bilinear form and the corresponding linear system,
#       applying any necessary transformations such as: eliminating boundary
#       conditions, applying conforming constraints for non-conforming AMR,
#       static condensation, etc.


if (static_cond):  a.EnableStaticCondensation()
a.Assemble();

A = mfem.SparseMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
## Here, original version calls hegith, which is not
## defined in the header...!?
print("Size of linear system: " + str(A.Size())) 

# 10. Solve 
M = mfem.GSSmoother(A)
mfem.PCG(A, M, B, X, 1, 500, 1e-12, 0.0);

# 11. Recover the solution as a finite element grid function.
a.RecoverFEMSolution(X, b, x);

# 12. Compute and print the L^2 norm of the error.
print("")
print("|| E_h - E ||_{L^2} = " + str(x.ComputeL2Error(E)))

mesh.PrintToFile('refined.mesh', 8)
x.SaveToFile('sol.gf', 8)



