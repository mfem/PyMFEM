'''
   MFEM example 12p

   See c++ version in the MFEM library for more detail 

   How to run:
      mpirun -np 2 python <arguments>
  
   Example of arguments:
      ex12p.py -m beam-tri.mesh
      ex12p.py -m beam-quad.mesh
      ex12p.py -m beam-tet.mesh -n 10 -o 2 -elast
      ex12p.py -m beam-hex.mesh
      ex12p.py -m beam-tri.mesh -o 2 -sys
      ex12p.py -m beam-quad.mesh -n 6 -o 3 -elast
      ex12p.py -m beam-quad-nurbs.mesh
      ex12p.py -m beam-hex-nurbs.mesh
'''
import sys
from os.path import expanduser, join, dirname
import numpy as np

from mfem.common.arg_parser import ArgParser
import mfem.par as mfem
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)

parser = ArgParser(description='Ex12 (linear elastisity)')
parser.add_argument('-m', '--mesh',
                    default='beam-tri.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument("-n", "--num-eigs",
                    action='store', default=5, type=int,
                    help="Number of desired eigenmodes.")
parser.add_argument("-s", "--seed",
                    action='store', default=66, type=int,
                    help="Random seed used to initialize LOBPCG.")
parser.add_argument('-elast', '--amg-for-elasticity',
                    action='store_true',
                    help='Use the special AMG elasticity solver (GM/LN approaches)',
                    dest='amg_elast', default=False)
parser.add_argument('-sys', '--amg-for-systems',
                    action='store_false',
                    help='Use  standard AMG for systems (unknown approach).',
                    dest='amg_elast', default=True)
parser.add_argument('-vis', '--visualization',
                    action='store_true', default=True,
                    help='Enable GLVis visualization')
args = parser.parse_args()

order = args.order
amg_elast = args.amg_elast
nev = args.num_eigs
visualization = args.visualization
seed = args.seed
if (myid == 0):
    parser.print_options(args)

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

# 3. Read the mesh from the given mesh file on all processors. We can handle
#    triangular, quadrilateral, tetrahedral, hexahedral, surface and volume
#    meshes with the same code
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
if (mesh.attributes.Max() < 2):
    if (myid == 0):
        print("Input mesh should have at least two materials!")
        print(" (See schematic in ex12p.cpp)")
    sys.exit()

# 4. Select the order of the finite element discretization space. For NURBS
#    meshes, we increase the order by degree elevation.
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

# 7. Define a parallel finite element space on the parallel mesh. Here we
#    use vector finite elements, i.e. dim copies of a scalar finite element
#    space. We use the ordering by vector dimension (the last argument of
#    the FiniteElementSpace constructor) which is expected in the systems
#    version of BoomerAMG preconditioner. For NURBS meshes, we use the
#    (degree elevated) NURBS space associated with the mesh nodes.
use_nodal_fespace = pmesh.NURBSext and not amg_elast
if use_nodal_fespace:
    fespace = pmesh.GetNodes().FESpace()
else:
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec, dim,
                                         mfem.Ordering.byVDIM)
size = fespace.GlobalTrueVSize()
if (myid == 0):
    print("Number of unknowns: " + str(size))
    print("Assembling")

# 8. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
#    element space corresponding to the linear elasticity integrator with
#    piece-wise constants coefficient lambda and mu, a simple mass matrix
#    needed on the right hand side of the generalized eigenvalue problem
#    below. The boundary conditions are implemented by marking only boundary
#    attribute 1 as essential. We use special values on the diagonal to
#    shift the Dirichlet eigenvalues out of the computational range. After
#    serial/parallel assembly we extract the corresponding parallel matrices
#    A and M.
lamb = mfem.Vector(pmesh.attributes.Max())
lamb.Assign(1.0)
lamb[0] = lamb[1]*50
lambda_func = mfem.PWConstCoefficient(lamb)
mu = mfem.Vector(pmesh.attributes.Max())
mu.Assign(1.0)
mu[0] = mu[1]*50
mu_func = mfem.PWConstCoefficient(mu)

ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
ess_bdr.Assign(0)
ess_bdr[0] = 1

a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_func, mu_func))
if (myid == 0):
    print("matrix ... ")
a.Assemble()
a.EliminateEssentialBCDiag(ess_bdr, 1.0)
a.Finalize()

m = mfem.ParBilinearForm(fespace)
m.AddDomainIntegrator(mfem.VectorMassIntegrator())
m.Assemble()

# shift the eigenvalue corresponding to eliminated dofs to a large value
m.EliminateEssentialBCDiag(ess_bdr,  sys.float_info.min)
m.Finalize()
if (myid == 0):
    print("done. ")

A = a.ParallelAssemble()
M = m.ParallelAssemble()

# 9. Define and configure the LOBPCG eigensolver and the BoomerAMG
#    preconditioner for A to be used within the solver. Set the matrices
#    which define the generalized eigenproblem A x = lambda M x.
amg = mfem.HypreBoomerAMG(A)
amg.SetPrintLevel(0)
if (amg_elast):
    amg.SetElasticityOptions(fespace)
else:
    amg.SetSystemsOptions(dim)

lobpcg = mfem.HypreLOBPCG(MPI.COMM_WORLD)
lobpcg.SetNumModes(nev)
lobpcg.SetRandomSeed(seed)
lobpcg.SetPreconditioner(amg)
lobpcg.SetMaxIter(100)
lobpcg.SetTol(1e-8)
lobpcg.SetPrecondUsageMode(1)
lobpcg.SetPrintLevel(1)
lobpcg.SetMassMatrix(M)
lobpcg.SetOperator(A)

# 10. Compute the eigenmodes and extract the array of eigenvalues. Define a
#     parallel grid function to represent each of the eigenmodes returned by
#     the solver.
eigenvalues = mfem.doubleArray()
lobpcg.Solve()
lobpcg.GetEigenvalues(eigenvalues)
x = mfem.ParGridFunction(fespace)


# 11. For non-NURBS meshes, make the mesh curved based on the finite element
#     space. This means that we define the mesh elements through a fespace
#     based transformation of the reference element. This allows us to save
#     the displaced mesh as a curved mesh when using high-order finite
#     element displacement field. We assume that the initial mesh (read from
#     the file) is not higher order curved mesh compared to the chosen FE
#     space.
if not use_nodal_fespace:
    pmesh.SetNodalFESpace(fespace)

# 12. Save the refined mesh and the modes in parallel. This output can be
#     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
pmesh.Print(mesh_name, 8)

for i in range(nev):
    x.Assign(lobpcg.GetEigenvector(i))
    sol_name = "mode_"+str(i).zfill(2)+"."+smyid
    x.Save(sol_name, 8)

# 13. Send the above data by socket to a GLVis server. Use the "n" and "b"
#     keys in GLVis to visualize the displacements.
if (visualization):
    mode_sock = mfem.socketstream("localhost", 19916)
    mode_sock.precision(8)

    for i in range(nev):
        if (myid == 0):
            print("Eigenmode " + str(i+1) + '/' + str(nev) +
                  ", Lambda = " + str(eigenvalues[i]))

        # convert eigenvector from HypreParVector to ParGridFunction
        x.Assign(lobpcg.GetEigenvector(i))

        mode_sock.send_text("parallel " + str(num_procs) + " " + str(myid))
        mode_sock.send_solution(pmesh,   x)
        mode_sock.send_text("window_title 'Eigenmode " + str(i+1) + '/' +
                            str(nev) + ", Lambda = " + str(eigenvalues[i]) + "'")

        c = None
        if (myid == 0):
            from builtins import input
            c = input("press (q)uit or (c)ontinue --> ")
        c = MPI.COMM_WORLD.bcast(c, root=0)
        if (c != 'c'):
            break

    mode_sock.close()
