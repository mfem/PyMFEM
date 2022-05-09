'''
   MFEM example 14p 

   How to run:
      mpirun -np 4 python <arguments>

   Example of arguments:
      ex14p.py -m inline-quad.mesh -o 0
      ex14p.py -m star.mesh -o 2
      ex14p.py -m escher.mesh -s 1
      ex14p.py -m fichera.mesh -s 1 -k 1
      ex14p.py -m square-disc-p3.mesh -o 3
      ex14p.py -m square-disc-nurbs.mesh -o 1
      ex14p.py -m disc-nurbs.mesh -rs 4 -o 2 -s 1 -k 0
      ex14p.py -m pipe-nurbs.mesh -o 1
      ex14p.py -m inline-segment.mesh -rs 5
      ex14p.py -m amr-quad.mesh -rs 3
      ex14p.py -m amr-hex.mesh
'''
import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np

import mfem.par as mfem
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

parser = ArgParser(description='Ex14p')
parser.add_argument('-m', '--mesh',
                    default='star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-rs', '--refine-serial',
                    action='store', default=-1, type=int,
                    help="Number of times to refine the mesh uniformly before parallel")
parser.add_argument('-rp', '--refine-parallel',
                    action='store', default=2, type=int,
                    help="Number of times to refine the mesh uniformly after parallel")
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-s', '--sigma',
                    action='store', default=-1.0, type=float,
                    help='\n'.join(["One of the two DG penalty parameters, typically +1/-1."
                                    " See the documentation of class DGDiffusionIntegrator."]))
parser.add_argument('-k', '--kappa',
                    action='store', default=-1.0, type=float,
                    help='\n'.join(["One of the two DG penalty parameters, should be positve."
                                    " Negative values are replaced with (order+1)^2."]))
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')
args = parser.parse_args()
ser_ref_levels = args.refine_serial
par_ref_levels = args.refine_parallel
order = args.order
sigma = args.sigma
kappa = args.kappa
visualization = args.visualization
if (kappa < 0):
    kappa = (order+1.)*(order+1.)
    args.kappa = kappa
if (myid == 0):
    parser.print_options(args)

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

# 3. Read the (serial) mesh from the given mesh file on all processors. We
#    can handle triangular, quadrilateral, tetrahedral and hexahedral meshes
#    with the same code. NURBS meshes are projected to second order meshes.
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()

# 4. Refine the serial mesh on all processors to increase the resolution. In
#    this example we do 'ser_ref_levels' of uniform refinement. By default,
#    or if ser_ref_levels < 0, we choose it to be the largest number that
#    gives a final mesh with no more than 10,000 elements.

if ser_ref_levels < 0:
    ser_ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ser_ref_levels):
    mesh.UniformRefinement()

if (mesh.NURBSext):
    mesh.SetCurvature(max(order, 1))

# 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
#    this mesh further in parallel to increase the resolution. Once the
#    parallel mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
for x in range(par_ref_levels):
    pmesh.UniformRefinement()

# 6. Define a parallel finite element space on the parallel mesh. Here we
#    use discontinuous finite elements of the specified order >= 0.
fec = mfem.DG_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)
glob_size = fespace.GlobalTrueVSize()
if (myid == 0):
    print('Number of unknowns: ' + str(glob_size))

# 7. Set up the parallel linear form b(.) which corresponds to the
#    right-hand side of the FEM linear system.
b = mfem.ParLinearForm(fespace)
one = mfem.ConstantCoefficient(1.0)
zero = mfem.ConstantCoefficient(0.0)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
b.AddBdrFaceIntegrator(
    mfem.DGDirichletLFIntegrator(zero, one, sigma, kappa))
b.Assemble()

# 8. Define the solution vector x as a parallel finite element grid function
#    corresponding to fespace. Initialize x with initial guess of zero.
x = mfem.ParGridFunction(fespace)
x.Assign(0.0)

# 9. Set up the bilinear form a(.,.) on the finite element space
#    corresponding to the Laplacian operator -Delta, by adding the Diffusion
#    domain integrator and the interior and boundary DG face integrators.
#    Note that boundary conditions are imposed weakly in the form, so there
#    is no need for dof elimination. After serial and parallel assembly we
#    extract the corresponding parallel matrix A.
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
a.AddInteriorFaceIntegrator(mfem.DGDiffusionIntegrator(one, sigma, kappa))
a.AddBdrFaceIntegrator(mfem.DGDiffusionIntegrator(one, sigma, kappa))
a.Assemble()
a.Finalize()

# 10. Define the parallel (hypre) matrix and vectors representing a(.,.),
#     b(.) and the finite element approximation.
A = a.ParallelAssemble()
B = b.ParallelAssemble()
X = x.ParallelProject()  # HypreParVector

del a
del b

# 11. Depending on the symmetry of A, define and apply a parallel PCG or
#     GMRES solver for AX=B using the BoomerAMG preconditioner from hypre.
amg = mfem.HypreBoomerAMG(A)
if sigma == -1.0:
    pcg = mfem.HyprePCG(A)
    pcg.SetTol(1e-12)
    pcg.SetMaxIter(200)
    pcg.SetPrintLevel(2)
    pcg.SetPreconditioner(amg)
    pcg.Mult(B, X)
else:
    gmres = mfem.GMRESSolver(MPI.COMM_WORLD)
    gmres.SetAbsTol(0.0)
    gmres.SetRelTol(1e-12)
    gmres.SetMaxIter(200)
    gmres.SetKDim(10)
    gmres.SetPrintLevel(1)
    gmres.SetOperator(A)
    gmres.SetPreconditioner(amg)
    gmres.Mult(B, X)
del amg

# 12. Extract the parallel grid function corresponding to the finite element
#     approximation X. This is the local solution on each processor.
x.Assign(X)

# 13. Save the refined mesh and the solution in parallel. This output can
#     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
sol_name = "sol."+smyid

pmesh.Print(mesh_name, 8)
x.Save(sol_name, 8)

# 14. Send the solution by socket to a GLVis server.
if (visualization):
    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock.send_text("parallel " + str(num_procs) + " " + str(myid))
    sol_sock.precision(8)
    sol_sock.send_solution(pmesh,  x)
