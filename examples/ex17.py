'''
   MFEM example 17 

   How to run:
      mpirun -np 2 python2.7 <arguments>

   Example of arguments:
       ex17.py -m beam-tri.mesh
       ex17.py -m beam-quad.mesh
       ex17.py -m beam-tet.mesh
       ex17.py -m beam-hex.mesh
       ex17.py -m beam-quad.mesh -r 2 -o 3
       ex17.py -m beam-quad.mesh -r 2 -o 2 -a 1 -k 1
       ex17.py -m beam-hex.mesh -r 2 -o 2

'''
import sys
import argparse
from os.path import expanduser, join
import numpy as np
from mfem import path

import mfem.ser as mfem
from mfem.ser import intArray


parser = argparse.ArgumentParser(description='Ex17')
parser.add_argument('-m', '--mesh',
                    default = 'beam-tri.mesh',
                    action = 'store', type = str,
                    help='Mesh file to use.')
parser.add_argument('-r', '--refine',
                    default = 'beam-tri.mesh',
                    action = 'store', default = -1, type=int,
                    help = "Number of times to refine the mesh uniformly, -1 for auto.")
parser.add_argument('-o', '--order',
                    action = 'store', default = 1, type=int,
                    help = "Finite element order (polynomial degree)");
parser.add_argument('-a', '--alpha',
                    action = 'store', default = -1.0, type=float,
          help = '\n'.join(["One of the two DG penalty parameters, typically +1/-1."
                            " See the documentation of class DGElasticityIntegrator."]))
parser.add_argument('-k', '--kappa',
                    action = 'store', default = -1.0, type=float,
          help = '\n'.join(["One of the two DG penalty parameters, should be positve."
                            " Negative values are replaced with (order+1)^2."]))
parser.add_argument('-vis', '--visualization',
                    action = 'store_true',
                    help='Enable GLVis visualization')

args = parser.parse_args()
ref_levels = args.refine
order = args.order
alpha = args.alpha;
kappa = args.kappa;
visualization = args.visualization
if (kappa < 0): kappa = (order+1)*(order+1)

# 
class InitDisplacement(mfem.VectorPyCoefficient):
    def __init__(self, dim):
       self.dim = dim
       mfem.VectorPyCoefficient.__init__(self, dim)
    def EvalValue(self, x):
       u = [0.0]*dim
       u[-1]  = -0.2*x[0]
       return tuple(u)

# 2. Read the mesh from the given mesh file.
meshfile =expanduser(join(path, 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1,1)
dim = mesh.Dimension()
if (mesh.attributes.Max() < 2 or 
    mesh.bdr_attributes.Max() < 2):
    print("\n".join(["Input mesh should have at least two materials and ", "two boundary attributes! (See schematic in ex17.cpp)\n"]))
    sys.exit()

# 3. Refine the mesh to increase the resolution.
ref_levels = int(np.floor(np.log(5000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
   mesh.UniformRefinement();

# Since NURBS meshes do not support DG integrators, we convert them to
# regular polynomial mesh of the specified (solution) order.
if (mesh.NURBSext):  mesh.SetCurvature(order)

# 4. Define a DG vector finite element space on the mesh. Here, we use
#    Gauss-Lobatto nodal basis because it gives rise to a sparser matrix
#    compared to the default Gauss-Legendre nodal basis.
fec = mfem.DG_FECollectio(order, dim, mfem.BasisType.GaussLobatto)
fespace = mfem.FiniteElementSpace(mesh, fec, dim)
print('Number of finite element unknowns: '+ str(fespace.GetVSize()))
print('Assembling:')

# 5. In this example, the Dirichlet boundary conditions are defined by
#    marking boundary attributes 1 and 2 in the marker Array 'dir_bdr'.
#    These b.c. are imposed weakly, by adding the appropriate boundary
#    integrators over the marked 'dir_bdr' to the bilinear and linear forms.
#    With this DG formulation, there are no essential boundary conditions.
ess_tdof_list = intArray()
dir_bdr = intArray(mesh.bdr_attributes.Max())
dir_bdr.Assign(0)
dir_bdr[0] = 1 # boundary attribute 1 is Dirichlet
dir_bdr[1] = 1 # boundary attribute 2 is Dirichlet

fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

# 6. Define the DG solution vector 'x' as a finite element grid function
#    corresponding to fespace. Initialize 'x' using the 'InitDisplacement'
#    function.
x = mfem.GridFunction(fespace) 
init_x = InitDisplacement(dim)
x.ProjectCoefficient(init_x)

# 7. Set up the Lame constants for the two materials. They are defined as
#    piece-wise (with respect to the element attributes) constant
#    coefficients, i.e. type PWConstCoefficient.
lamb = mfem.Vector(mesh.attributes.Max())  # lambda is not possible in python
lamb.Assign(1.0)
lamb[0] = 50.
lambda_c = mfem.PWConstCoefficient(lamb)
mu = mfem.Vector(mesh.attributes.Max())
mu.Assign(1.0);
mu[0] = 50.0
mu_c = mfem.PWConstCoefficient(mu)

# 8. Set up the linear form b(.) which corresponds to the right-hand side of
#    the FEM linear system. In this example, the linear form b(.) consists
#    only of the terms responsible for imposing weakly the Dirichlet
#    boundary conditions, over the attributes marked in 'dir_bdr'. The
#    values for the Dirichlet boundary condition are taken from the
#    VectorFunctionCoefficient 'x_init' which in turn is based on the
#    function 'InitDisplacement'.
b = mfem.LinearForm(fespace)
print('r.h.s ...')
integrator = mfem.DGElasticityDirichletLFIntegrator(init_x, lambda_c, mu_c, alpha, kappa)
b.AddBdrFaceIntegrator(integrator , dir_bdr)
b.Assemble()

# 9. Set up the bilinear form a(.,.) on the DG finite element space
#    corresponding to the linear elasticity integrator with coefficients
#    lambda and mu as defined above. The additional interior face integrator
#    ensures the weak continuity of the displacement field. The additional
#    boundary face integrator works together with the boundary integrator
#    added to the linear form b(.) to impose weakly the Dirichlet boundary
#    conditions.
a = mfem.BilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_c, mu_c))
a.AddInteriorFaceIntegrator(mfem.DGElasticityIntegrator(lambda_c, mu_c, alpha, kappa))
a.AddBdrFaceIntegrator(mfem.DGElasticityIntegrator(lambda_c, mu_c, alpha, kappa), dir_bdr)
print('matrix ...')
a.Assemble()

A = mfem.SparseMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
print('...done')

A.PrintInfo(sys.stdout)

# 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
#     solve the system Ax=b with PCG for the symmetric formulation, or GMRES
#     for the non-symmetric.
M = mfem.GSSmoother(A)
rtol = 1e-6
if (alpha == -1.0):
    mfem.PCG(A, M, B, X, 3, 5000, rtol*rtol, 0.0)
else:
    mfem.GMRES(A, M, B, X, 3, 5000, 50, rtol*rtol, 0.0)

# 12. Recover the solution as a finite element grid function 'x'.
a.RecoverFEMSolution(X, b, x)

# 13. Use the DG solution space as the mesh nodal space. This allows us to
#     save the displaced mesh as a curved DG mesh.
mesh.SetNodalFESpace(fespace)
if (visualization):
     reference_nodes = mesh.GetNodes()
# 14. Save the displaced mesh and minus the solution (which gives the
#     backward displacements to the reference mesh). This output can be
#     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
nodes = mesh.GetNodes()
nodes += x
x.Neg()
mesh.PrintToFile('displaced.mesh', 8)
x.SaveToFile('sol.gf', 8)

# 15. Visualization: send data by socket to a GLVis server.
