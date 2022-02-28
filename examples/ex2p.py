'''
   MFEM example 2

   See c++ version in the MFEM library for more detail 

   How to run:
      mpirun -np 2 python2.7 <arguments>
  
   Example of arguments:
      ex2p.py -m star.mesh
      ex2p.py -m square-disc.mesh
      ex2p.py -m escher.mesh
      ex2p.py -m fichera.mesh
      ex2p.py -m beam-tri.mesh -o 2 -sys   
      ex2p.py -m beam-quad.mesh -o 3 -elast
      ex2p.py -m beam-quad.mesh -o 3 -sc
'''
import numpy as np
from os.path import expanduser, join, dirname
import sys
from mfem import path
from mfem.common.arg_parser import ArgParser
import mfem.par as mfem
from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


parser = ArgParser(description='Ex2 (linear elastisity)')
parser.add_argument('-m', '--mesh',
                    default='beam-tri.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree) or -1 for isoparametric space.")
parser.add_argument('-sc', '--static-condensation',
                    action='store_false',
                    help="Enable static condensation.")
parser.add_argument('-elast', '--amg-for-elasticity',
                    action='store_true',
                    help='Use the special AMG elasticity solver (GM/LN approaches)',
                    dest='amg_elast', default=False)
parser.add_argument('-sys', '--amg-for-systems',
                    action='store_false',
                    help='Use  standard AMG for systems (unknown approach).',
                    dest='amg_elast', default=True)
args = parser.parse_args()

order = args.order
static_cond = args.static_condensation
amg_elast = args.amg_elast
visualization = args.visualization
if (myid == 0):
    parser.print_options(args)

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

#  3. Read the (serial) mesh from the given mesh file on all processors.  We
#     can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
#     and volume meshes with the same code.
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
if (mesh.attributes.Max() < 2 or mesh.bdr_attributes.Max() < 2):
    if (myid == 0):
        print('\n'.join(['Input mesh should have at least two materials and',
                         'two boundary attributes! (See schematic in ex2.cpp)']))
    sys.exit()

#  4. Select the order of the finite element discretization space. For NURBS
#     meshes, we increase the order by degree elevation.
if (mesh.NURBSext and order > mesh.NURBSext.GetOrder()):
    mesh.DegreeElevate(order - mesh.NURBSext.GetOrder())

#  5. Refine the serial mesh on all processors to increase the resolution. In
#     this example we do 'ref_levels' of uniform refinement. We choose
#     'ref_levels' to be the largest number that gives a final mesh with no
#     more than 1,000 elements.
ref_levels = int(np.floor(np.log(1000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

#  6. Define a parallel mesh by a partitioning of the serial mesh. Refine
#     this mesh further in parallel to increase the resolution. Once the
#     parallel mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
par_ref_levels = 1
for l in range(par_ref_levels):
    pmesh.UniformRefinement()

#   7. Define a parallel finite element space on the parallel mesh. Here we
#      use vector finite elements, i.e. dim copies of a scalar finite element
#      space. We use the ordering by vector dimension (the last argument of
#      the FiniteElementSpace constructor) which is expected in the systems
#      version of BoomerAMG preconditioner. For NURBS meshes, we use the
#      (degree elevated) NURBS space associated with the mesh nodes.
use_nodal_fespace = pmesh.NURBSext and not amg_elast
if use_nodal_fespace:
    fespace = pmesh.GetNodes().FESpace()
else:
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec, dim,
                                         mfem.Ordering.byVDIM)

size = fespace.GlobalTrueVSize()
if (myid == 0):
    print("Number of finite element unknowns: " + str(size))

# 8. Determine the list of true (i.e. parallel conforming) essential
#    boundary dofs. In this example, the boundary conditions are defined by
#    marking only boundary attribute 1 from the mesh as essential and
#    converting it to a list of true dofs.
ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
ess_tdof_list = mfem.intArray()
ess_bdr.Assign(0)
ess_bdr[0] = 1
fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)
print("here")

#  9. Set up the parallel linear form b(.) which corresponds to the
#     right-hand side of the FEM linear system. In this case, b_i equals the
#     boundary integral of f*phi_i where f represents a "pull down" force on
#     the Neumann part of the boundary and phi_i are the basis functions in
#     the finite element fespace. The force is defined by the object f, which
#     is a vector of Coefficient objects. The fact that f is non-zero on
#     boundary attribute 2 is indicated by the use of piece-wise constants
#     coefficient for its last component.
f = mfem.VectorArrayCoefficient(dim)
for i in range(dim-1):
    f.Set(i, mfem.ConstantCoefficient(0.0))

pull_force = mfem.Vector([0]*pmesh.bdr_attributes.Max())
pull_force[1] = -1.0e-2
f.Set(dim-1, mfem.PWConstCoefficient(pull_force))


b = mfem.ParLinearForm(fespace)
b.AddBoundaryIntegrator(mfem.VectorBoundaryLFIntegrator(f))
print('r.h.s. ...')
b.Assemble()

# 10. Define the solution vector x as a parallel finite element grid
#     function corresponding to fespace. Initialize x with initial guess of
#     zero, which satisfies the boundary conditions.
x = mfem.GridFunction(fespace)
x.Assign(0.0)
print('here')
# 11. Set up the parallel bilinear form a(.,.) on the finite element space
#     corresponding to the linear elasticity integrator with piece-wise
#     constants coefficient lambda and mu.

lamb = mfem.Vector(pmesh.attributes.Max())
lamb.Assign(1.0)
lamb[0] = lamb[1]*50
lambda_func = mfem.PWConstCoefficient(lamb)
mu = mfem.Vector(pmesh.attributes.Max())
mu.Assign(1.0)
mu[0] = mu[1]*50
mu_func = mfem.PWConstCoefficient(mu)
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_func, mu_func))

#  12. Assemble the parallel bilinear form and the corresponding linear
#      system, applying any necessary transformations such as: parallel
#      assembly, eliminating boundary conditions, applying conforming
#      constraints for non-conforming AMR, static condensation, etc.
if (myid == 0):
    print('matrix...')
if (static_cond):
    a.EnableStaticCondensation()
a.Assemble()

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
if (myid == 0):
    print('...done')
    print("Size of linear system: " + str(A.GetGlobalNumRows()))

#  13. Define and apply a parallel PCG solver for A X = B with the BoomerAMG
#      preconditioner from hypre.
amg = mfem.HypreBoomerAMG(A)
if (amg_elast and not a.StaticCondensationIsEnabled()):
    amg.SetElasticityOptions(fespace)
else:
    amg.SetSystemsOptions(dim)
pcg = mfem.HyprePCG(A)
pcg.SetTol(1e-8)
pcg.SetMaxIter(500)
pcg.SetPrintLevel(2)
pcg.SetPreconditioner(amg)
pcg.Mult(B, X)

#  14. Recover the parallel grid function corresponding to X. This is the
#      local finite element solution on each processor.
a.RecoverFEMSolution(X, b, x)

#  15. For non-NURBS meshes, make the mesh curved based on the finite element
#      space. This means that we define the mesh elements through a fespace
#      based transformation of the reference element.  This allows us to save
#      the displaced mesh as a curved mesh when using high-order finite
#      element displacement field. We assume that the initial mesh (read from
#      the file) is not higher order curved mesh compared to the chosen FE
#      space.
if (not use_nodal_fespace):
    pmesh.SetNodalFESpace(fespace)

#  16. Save in parallel the displaced mesh and the inverted solution (which
#      gives the backward displacements to the original grid). This output
#      can be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".

nodes = pmesh.GetNodes()
nodes += x
x *= -1

smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
sol_name = "sol."+smyid

pmesh.Print(mesh_name, 8)
x.Save(sol_name, 8)

#  17. Send the above data by socket to a GLVis server.  Use the "n" and "b"
#      keys in GLVis to visualize the displacements.

if (visualization):
    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock.send_text("parallel " + str(num_procs) + " " + str(myid))
    sol_sock.precision(8)
    sol_sock.send_solution(pmesh,  x)
