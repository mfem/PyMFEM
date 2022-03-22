'''
   MFEM example 7

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from mpi4py import MPI
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, exp


num_proc = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

elem_type = 1
ref_levels = 2
par_ref_levels = 2
amr = 0
order = 2
always_snap = False

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

if elem_type == 1:
    Nvert = 8
    Nelem = 6
else:
    Nvert = 6
    Nelem = 8


mesh = mfem.Mesh(2, Nvert, Nelem, 0, 3)

if elem_type == 0:
    tri_v = [[1.,  0.,  0.], [0.,  1.,  0.], [-1.,  0.,  0.],
             [0., -1.,  0.], [0.,  0.,  1.], [0.,  0., -1.]]
    tri_e = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4],
             [1, 0, 5], [2, 1, 5], [3, 2, 5], [0, 3, 5]]
    for j in range(Nvert):
        mesh.AddVertex(tri_v[j])
    for j in range(Nelem):
        mesh.AddTriangle(tri_e[j], j+1)
    mesh.FinalizeTriMesh(1, 1, True)
else:
    quad_v = [[-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
              [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1]]

    quad_e = [[3, 2, 1, 0], [0, 1, 5, 4], [1, 2, 6, 5],
              [2, 3, 7, 6], [3, 0, 4, 7], [4, 5, 6, 7]]
    for j in range(Nvert):
        mesh.AddVertex(quad_v[j])
    for j in range(Nelem):
        mesh.AddQuad(quad_e[j], j+1)
    mesh.FinalizeQuadMesh(1, 1, True)

#  Set the space for the high-order mesh nodes.
fec = mfem.H1_FECollection(order, mesh.Dimension())
nodal_fes = mfem.FiniteElementSpace(mesh, fec, mesh.SpaceDimension())
mesh.SetNodalFESpace(nodal_fes)


def SnapNodes(mesh):
    nodes = mesh.GetNodes()
    node = mfem.Vector(mesh.SpaceDimension())

    for i in np.arange(nodes.FESpace().GetNDofs()):
        for d in np.arange(mesh.SpaceDimension()):
            node[d] = nodes[nodes.FESpace().DofToVDof(i, d)]

        node /= node.Norml2()
        for d in range(mesh.SpaceDimension()):
            nodes[nodes.FESpace().DofToVDof(i, d)] = node[d]


#  3. Refine the mesh while snapping nodes to the sphere.
for l in range(ref_levels+1):
    if l > 0:
        mesh.UniformRefinement()
    if always_snap:
        SnapNodes(mesh)
if amr == 1:
    for l in range(5):
        mesh.RefineAtVertex(mfem.Vertex(0, 0, 1))
    SnapNodes(mesh)
elif amr == 2:
    for l in range(4):
        mesh.RandomRefinement(0.5)
    SnapNodes(mesh)

pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
for i in range(par_ref_levels):
    pmesh.UniformRefinement()
    if always_snap:
        SnapNodes(pmesh)
if (not always_snap) or (par_ref_levels < 1):
    SnapNodes(pmesh)
if amr == 1:
    for l in range(5):
        pmesh.RefineAtVertex(mfem.Vertex(0, 0, 1))
    SnapNodes(pmesh)
elif amr == 2:
    for l in range(4):
        pmesh.RandomRefinement(0.5)
    SnapNodes(pmesh)


#  4. Define a finite element space on the mesh. Here we use isoparametric
#      finite elements -- the same as the mesh nodes.
fespace = mfem.ParFiniteElementSpace(pmesh, fec)
size = fespace.GlobalTrueVSize()
if myid == 0:
    print("Number of unknowns: " + str(size))

#   5. Set up the linear form b(.) which corresponds to the right-hand side of
#      the FEM linear system, which in this case is (1,phi_i) where phi_i are
#      the basis functions in the finite element fespace.


class analytic_rhs(mfem.PyCoefficient):
    def EvalValue(self, x):
        l2 = np.sum(x**2)
        return 7.*x[0]*x[1]/l2


class analytic_solution(mfem.PyCoefficient):
    def EvalValue(self, x):
        l2 = np.sum(x**2)
        return x[0]*x[1]/l2


b = mfem.ParLinearForm(fespace)
one = mfem.ConstantCoefficient(1.0)
rhs_coef = analytic_rhs()
sol_coef = analytic_solution()
b.AddDomainIntegrator(mfem.DomainLFIntegrator(rhs_coef))
b.Assemble()

#  6. Define the solution vector x as a finite element grid function
#     corresponding to fespace. Initialize x with initial guess of zero.
x = mfem.ParGridFunction(fespace)
x.Assign(0.0)

#  7. Set up the bilinear form a(.,.) on the finite element space
#     corresponding to the Laplacian operator -Delta, by adding the Diffusion
#     and Mass domain integrators.
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
a.AddDomainIntegrator(mfem.MassIntegrator(one))

#  8. Assemble the linear system, apply conforming constraints, etc.
a.Assemble()
A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
empty_tdof_list = mfem.intArray()

a.FormLinearSystem(empty_tdof_list, x, b, A, X, B)

amg = mfem.HypreBoomerAMG(A)
pcg = mfem.HyprePCG(A)
pcg.SetTol(1e-12)
pcg.SetMaxIter(200)
pcg.SetPrintLevel(2)
pcg.SetPreconditioner(amg)
pcg.Mult(B, X)

# 10. Recover the solution as a finite element grid function.
a.RecoverFEMSolution(X, b, x)

# 12. Compute and print the L^2 norm of the error.
err = x.ComputeL2Error(sol_coef)
if myid == 0:
    print("L2 norm of error = " + str(err))


smyid = '{:0>6d}'.format(myid)
mesh_name = "sphere_refined."+smyid
sol_name = "sol."+smyid

pmesh.Print(mesh_name, 8)
x.Save(sol_name, 8)
