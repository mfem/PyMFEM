'''
   MFEM example 11

   See c++ version in the MFEM library for more detail 

   How to run:
      mpirun -np 2 python <arguments>
  
   Example of arguments:
      ex11p.py -m square-disc.mesh
      ex11p.py -m star.mesh
      ex11p.py -m escher.mesh
      ex11p.py -m fichera.mesh
      ex11p.py -m square-disc-p2.vtk -o 2
      ex11p.py -m square-disc-p3.mesh -o 3
      ex11p.py -m square-disc-nurbs.mesh -o -1
      ex11p.py -m disc-nurbs.mesh -o -1 -n 20
      ex11p.py -m pipe-nurbs.mesh -o -1
      ex11p.py -m ball-nurbs.mesh -o 2
      ex11p.py -m star-surf.mesh
      ex11p.py -m square-disc-surf.mesh
      ex11p.py -m inline-segment.mesh
      ex11p.py -m amr-quad.mesh
      ex11p.py -m amr-hex.mesh
      ex11p.py -m mobius-strip.mesh -n 8
      ex11p.py -m klein-bottle.mesh -n 10
'''
import sys
from os.path import expanduser, join, dirname
import numpy as np

from mfem.common.arg_parser import ArgParser
import mfem.par as mfem
from mpi4py import MPI


num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


parser = ArgParser(description='Ex11 ')
parser.add_argument('-m', '--mesh',
                    default='star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-rs', '--refine-serial',
                    default=2,
                    action='store', type=int,
                    help="Number of times to refine the mesh uniformly in serial.")

parser.add_argument('-rp', '--refine-parallel',
                    default=1,
                    action='store', type=int,
                    help="Number of times to refine the mesh uniformly in parallel.")
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help=("Finite element order (polynomial degree) or -1 for isoparametric space."))
parser.add_argument("-n", "--num-eigs",
                    action='store', default=5, type=int,
                    help="Number of desired eigenmodes.")
parser.add_argument("-sp", "--strumpack",
                    action='store_true', default=False,
                    help="Use the STRUMPACK Solver.")
parser.add_argument('-vis', '--visualization',
                    action='store_true', default=True,
                    help='Enable GLVis visualization')
args = parser.parse_args()


ser_ref_levels = args.refine_serial
par_ref_levels = args.refine_parallel
order = args.order
nev = args.num_eigs
visualization = args.visualization
use_strumpack = args.strumpack
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

# 4. Refine the serial mesh on all processors to increase the resolution. In
#    this example we do 'ref_levels' of uniform refinement (2 by default, or
#    specified on the command line with -rs).
for x in range(ser_ref_levels):
    mesh.UniformRefinement()

# 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
#    this mesh further in parallel to increase the resolution (1 time by
#    default, or specified on the command line with -rp). Once the parallel
#    mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
for l in range(par_ref_levels):
    pmesh.UniformRefinement()

# 6. Define a parallel finite element space on the parallel mesh. Here we
#    use continuous Lagrange finite elements of the specified order. If
#    order < 1, we instead use an isoparametric/isogeometric space.
if order > 0:
    fec = mfem.H1_FECollection(order, dim)
elif pmesh.GetNodes():
    fec = pmesh.GetNodes().OwnFEC()
else:
    fec = mfem.H1_FECollection(1, dim)

fespace = mfem.ParFiniteElementSpace(pmesh, fec)
fe_size = fespace.GlobalTrueVSize()

if (myid == 0):
    print('Number of unknowns: ' + str(fe_size))

# 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
#    element space. The first corresponds to the Laplacian operator -Delta,
#    while the second is a simple mass matrix needed on the right hand side
#    of the generalized eigenvalue problem below. The boundary conditions
#    are implemented by elimination with special values on the diagonal to
#    shift the Dirichlet eigenvalues out of the computational range. After
#    serial and parallel assembly we extract the corresponding parallel
#    matrices A and M.
one = mfem.ConstantCoefficient(1.0)

ess_bdr = mfem.intArray()
if pmesh.bdr_attributes.Size() != 0:
    ess_bdr.SetSize(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)

a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
if pmesh.bdr_attributes.Size() == 0:
    # Add a mass term if the mesh has no boundary, e.g. periodic mesh or
    # closed surface.
    a.AddDomainIntegrator(mfem.MassIntegrator(one))

a.Assemble()
a.EliminateEssentialBCDiag(ess_bdr, 1.0)
a.Finalize()


m = mfem.ParBilinearForm(fespace)
m.AddDomainIntegrator(mfem.MassIntegrator(one))
m.Assemble()

# shift the eigenvalue corresponding to eliminated dofs to a large value
m.EliminateEssentialBCDiag(ess_bdr,  3.0e-300)
m.Finalize()

A = a.ParallelAssemble()
M = m.ParallelAssemble()

if use_strumpack:
    import mfem.par.strumpack as strmpk
    Arow = strmpk.STRUMPACKRowLocMatrix(A)


# 8. Define and configure the LOBPCG eigensolver and the BoomerAMG
#    preconditioner for A to be used within the solver. Set the matrices
#    which define the generalized eigenproblem A x = lambda M x.
#    We don't support SuperLU

if use_strumpack:
    args = ["--sp_hss_min_sep_size", "128", "--sp_enable_hss"]
    strumpack = strmpk.STRUMPACKSolver(args, MPI.COMM_WORLD)
    strumpack.SetPrintFactorStatistics(True)
    strumpack.SetPrintSolveStatistics(False)
    strumpack.SetKrylovSolver(strmpk.KrylovSolver_DIRECT)
    strumpack.SetReorderingStrategy(strmpk.ReorderingStrategy_METIS)
    strumpack.SetMC64Job(strmpk.MC64Job_NONE)
    # strumpack.SetSymmetricPattern(True)
    strumpack.SetOperator(Arow)
    strumpack.SetFromCommandLine()
    precond = strumpack
else:
    amg = mfem.HypreBoomerAMG(A)
    amg.SetPrintLevel(0)
    precond = amg

lobpcg = mfem.HypreLOBPCG(MPI.COMM_WORLD)
lobpcg.SetNumModes(nev)
lobpcg.SetPreconditioner(precond)
lobpcg.SetMaxIter(200)
lobpcg.SetTol(1e-8)
lobpcg.SetPrecondUsageMode(1)
lobpcg.SetPrintLevel(1)
lobpcg.SetMassMatrix(M)
lobpcg.SetOperator(A)

# 9. Compute the eigenmodes and extract the array of eigenvalues. Define a
#    parallel grid function to represent each of the eigenmodes returned by
#    the solver.
eigenvalues = mfem.doubleArray()
lobpcg.Solve()
lobpcg.GetEigenvalues(eigenvalues)
x = mfem.ParGridFunction(fespace)

# 10. Save the refined mesh and the modes in parallel. This output can be
#     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".

smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
pmesh.Print(mesh_name, 8)

for i in range(nev):
    x.Assign(lobpcg.GetEigenvector(i))
    sol_name = "mode_"+str(i).zfill(2)+"."+smyid
    x.Save(sol_name, 8)

# 11. Send the solution by socket to a GLVis server.
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
