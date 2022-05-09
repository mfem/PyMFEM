'''
   MFEM example 8p 

   See c++ version in the MFEM library for more detail 

   How to run:
      mpirun -np 2 python <arguments>

   Example of arguments:
      ex8p.py -m square-disc.mesh
      ex8p.py -m star.mesh
      ex8p.py -m escher.mesh
      ex8p.py -m fichera.mesh
      ex8p.py -m square-disc-p3.mesh
      ex8p.py -m star-surf.mesh -o 2
'''
import sys
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt

from mfem.common.arg_parser import ArgParser

import mfem.par as mfem
from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


parser = ArgParser(description='Ex8p')
parser.add_argument('-m', '--mesh',
                    default='star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-vis', '--visualization',
                    action='store_true', default=True,
                    help='Enable GLVis visualization')

args = parser.parse_args()

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

order = args.order
visualization = args.visualization
if myid == 0:
    parser.print_options(args)

# 3. Read the (serial) mesh from the given mesh file on all processors.  We
#    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
#    and volume meshes with the same code.

meshfile = expanduser(join(dirname(__file__), '..', 'data', 'star.mesh'))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()

# 4. Refine the serial mesh on all processors to increase the resolution. In
#    this example we do 'ref_levels' of uniform refinement. We choose
#    'ref_levels' to be the largest number that gives a final mesh with no
#    more than 10,000 elements.
ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

# 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
#    this mesh further in parallel to increase the resolution. Once the
#    parallel mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
par_ref_levels = 1
for l in range(par_ref_levels):
    pmesh.UniformRefinement()
pmesh.ReorientTetMesh()

# 6. Define the trial, interfacial (trace) and test DPG spaces:
#    - The trial space, x0_space, contains the non-interfacial unknowns and
#      has the essential BC.
#    - The interfacial space, xhat_space, contains the interfacial unknowns
#      and does not have essential BC.
#    - The test space, test_space, is an enriched space where the enrichment
#      degree may depend on the spatial dimension of the domain, the type of
#      the mesh and the trial space order.
trial_order = order
trace_order = order - 1
test_order = order  # reduced order, full order is (order + dim - 1)

if (dim == 2 and (order % 2 == 0 or (pmesh.MeshGenerator() & 2 and order > 1))):
    test_order = test_order + 1
if (test_order < trial_order):
    if myid == 0:
        print("Warning, test space not enriched enough to handle primal trial space")

x0_fec = mfem.H1_FECollection(trial_order, dim)
xhat_fec = mfem.RT_Trace_FECollection(trace_order, dim)
test_fec = mfem.L2_FECollection(test_order, dim)

x0_space = mfem.ParFiniteElementSpace(pmesh, x0_fec)
xhat_space = mfem.ParFiniteElementSpace(pmesh, xhat_fec)
test_space = mfem.ParFiniteElementSpace(pmesh, test_fec)

glob_true_s0 = x0_space.GlobalTrueVSize()
glob_true_s1 = xhat_space.GlobalTrueVSize()
glob_true_s_test = test_space.GlobalTrueVSize()

if myid == 0:
    print('\n'.join(["nNumber of Unknowns",
                     " Trial space,     X0   : " + str(glob_true_s0) +
                     " (order " + str(trial_order) + ")",
                     " Interface space, Xhat : " + str(glob_true_s1) +
                     " (order " + str(trace_order) + ")",
                     " Test space,      Y    : " + str(glob_true_s_test) +
                     " (order " + str(test_order) + ")"]))

# 7. Set up the linear form F(.) which corresponds to the right-hand side of
#    the FEM linear system, which in this case is (f,phi_i) where f=1.0 and
#    phi_i are the basis functions in the test finite element fespace.
one = mfem.ConstantCoefficient(1.0)
F = mfem.ParLinearForm(test_space)
F.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
F.Assemble()

x0 = mfem.ParGridFunction(x0_space)
x0.Assign(0.0)

# 8. Set up the mixed bilinear form for the primal trial unknowns, B0,
#    the mixed bilinear form for the interfacial unknowns, Bhat,
#    the inverse stiffness matrix on the discontinuous test space, Sinv,
#    and the stiffness matrix on the continuous trial space, S0.
ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
ess_bdr.Assign(1)
ess_dof = mfem.intArray()
x0_space.GetEssentialVDofs(ess_bdr, ess_dof)

B0 = mfem.ParMixedBilinearForm(x0_space, test_space)
B0.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
B0.Assemble()
B0.EliminateEssentialBCFromTrialDofs(ess_dof, x0, F)
B0.Finalize()

Bhat = mfem.ParMixedBilinearForm(xhat_space, test_space)
Bhat.AddTraceFaceIntegrator(mfem.TraceJumpIntegrator())
Bhat.Assemble()
Bhat.Finalize()

Sinv = mfem.ParBilinearForm(test_space)
Sum = mfem.SumIntegrator()
Sum.AddIntegrator(mfem.DiffusionIntegrator(one))
Sum.AddIntegrator(mfem.MassIntegrator(one))
Sinv.AddDomainIntegrator(mfem.InverseIntegrator(Sum))
Sinv.Assemble()
Sinv.Finalize()

S0 = mfem.ParBilinearForm(x0_space)
S0.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
S0.Assemble()
S0.EliminateEssentialBC(ess_bdr)
S0.Finalize()

matB0 = B0.ParallelAssemble()
del B0
matBhat = Bhat.ParallelAssemble()
del Bhat
matSinv = Sinv.ParallelAssemble()
del Sinv
matS0 = S0.ParallelAssemble()
del S0

# 9. Define the block structure of the problem, by creating the offset
#    variables. Also allocate two BlockVector objects to store the solution
#    and rhs.
x0_var = 0
xhat_var = 1
NVAR = 2  # enum in C

true_s0 = x0_space.TrueVSize()
true_s1 = xhat_space.TrueVSize()
true_s_test = test_space.TrueVSize()


true_offsets = mfem.intArray([0, true_s0, true_s0+true_s1])
true_offsets_test = mfem.intArray([0, true_s_test])


x = mfem.BlockVector(true_offsets)
b = mfem.BlockVector(true_offsets)
x.Assign(0.0)
b.Assign(0.0)

# 10. Set up the 1x2 block Least Squares DPG operator, B = [B0 Bhat],
#     the normal equation operator, A = B^t Sinv B, and
#     the normal equation right-hand-size, b = B^t Sinv F.

B = mfem.BlockOperator(true_offsets_test, true_offsets)
B.SetBlock(0, 0, matB0)
B.SetBlock(0, 1, matBhat)

A = mfem.RAPOperator(B, matSinv, B)

trueF = F.ParallelAssemble()

SinvF = mfem.HypreParVector(test_space)
matSinv.Mult(trueF, SinvF)
B.MultTranspose(SinvF, b)

# 11. Set up a block-diagonal preconditioner for the 2x2 normal equation
#
#        [ S0^{-1}     0     ]
#        [   0     Shat^{-1} ]      Shat = (Bhat^T Sinv Bhat)
#
#     corresponding to the primal (x0) and interfacial (xhat) unknowns.
#     Since the Shat operator is equivalent to an H(div) matrix reduced to
#     the interfacial skeleton, we approximate its inverse with one V-cycle
#     of the ADS preconditioner from the hypre library (in 2D we use AMS for
#     the rotated H(curl) problem).
S0inv = mfem.HypreBoomerAMG(matS0)
S0inv.SetPrintLevel(0)

Shat = mfem.RAP(matSinv, matBhat)
if (dim == 2):
    Shatinv = mfem.HypreAMS(Shat, xhat_space)
else:
    Shatinv = mfem.HypreADS(Shat, xhat_space)

P = mfem.BlockDiagonalPreconditioner(true_offsets)
P.SetDiagonalBlock(0, S0inv)
P.SetDiagonalBlock(1, Shatinv)

# 12. Solve the normal equation system using the PCG iterative solver.
#     Check the weighted norm of residual for the DPG least square problem.
#     Wrap the primal variable in a GridFunction for visualization purposes.
pcg = mfem.CGSolver(MPI.COMM_WORLD)
pcg.SetOperator(A)
pcg.SetPreconditioner(P)
pcg.SetRelTol(1e-6)
pcg.SetMaxIter(200)
pcg.SetPrintLevel(1)
pcg.Mult(b, x)

LSres = mfem.HypreParVector(test_space)
tmp = mfem.HypreParVector(test_space)

B.Mult(x, LSres)
LSres -= trueF
matSinv.Mult(LSres, tmp)

res = sqrt(mfem.InnerProduct(LSres, tmp))

if (myid == 0):
    print("\n|| B0*x0 + Bhat*xhat - F ||_{S^-1} = " + str(res))


x0.Distribute(x.GetBlock(x0_var))

# 13. Save the refined mesh and the solution in parallel. This output can
#     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
sol_name = "sol."+smyid
pmesh.Print(mesh_name, 8)
x0.Save(sol_name, 8)

# 14. Send the solution by socket to a GLVis server.
if visualization:
    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock.send_text("parallel " + str(num_procs) + " " + str(myid))
    sol_sock.precision(8)
    sol_sock.send_solution(pmesh, x0)
