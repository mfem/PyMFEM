'''
   MFEM example 8 (DPG method )

   See c++ version in the MFEM library for more detail 

   How to run:
      python <arguments>

   Example of arguments:
      ex8.py -m square-disc.mesh
      ex8.py -m star.mesh
      ex8.py -m escher.mesh
      ex8.py -m fichera.mesh
      ex8.py -m square-disc-p3.mesh
      ex8.py -m star-surf.mesh -o 2
      ex8.py -m mobius-strip.mesh
'''
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt

from mfem.common.arg_parser import ArgParser
from mfem import path
import mfem.ser as mfem

parser = ArgParser(description='Ex8')
parser.add_argument('-m', '--mesh',
                    default='../data/star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-vis', '--visualization',
                    action='store_true', default=True,
                    help='Enable GLVis visualization')

args = parser.parse_args()
order = args.order
visualization = args.visualization

#   2. Read the mesh from the given mesh file. We can handle triangular,
#      quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
#      the same code.
path = dirname((__file__))
meshfile = expanduser(join(path, '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)

dim = mesh.Dimension()

#   3. Refine the mesh to increase the resolution. In this example we do
#      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
#      largest number that gives a final mesh with no more than 10,000
#      elements.
ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

#  4. Define the trial, interfacial (trace) and test DPG spaces:
#     - The trial space, x0_space, contains the non-interfacial unknowns and
#       has the essential BC.
#     - The interfacial space, xhat_space, contains the interfacial unknowns
#       and does not have essential BC.
#     - The test space, test_space, is an enriched space where the enrichment
#       degree may depend on the spatial dimension of the domain, the type of
#       the mesh and the trial space order.

trial_order = order
trace_order = order - 1
test_order = order  # reduced order, full order is (order + dim - 1)

if (dim == 2 and (order % 2 == 0 or (mesh.MeshGenerator() & 2 and order > 1))):
    test_order = test_order + 1
if (test_order < trial_order):
    print("Warning, test space not enriched enough to handle primal trial space")

x0_fec = mfem.H1_FECollection(trial_order, dim)
xhat_fec = mfem.RT_Trace_FECollection(trace_order, dim)
test_fec = mfem.L2_FECollection(test_order, dim)

x0_space = mfem.FiniteElementSpace(mesh, x0_fec)
xhat_space = mfem.FiniteElementSpace(mesh, xhat_fec)
test_space = mfem.FiniteElementSpace(mesh, test_fec)

# 5. Define the block structure of the problem, by creating the offset
#    variables. Also allocate two BlockVector objects to store the solution
#    and rhs.

x0_var = 0
xhat_var = 1
NVAR = 2

s0 = x0_space.GetVSize()
s1 = xhat_space.GetVSize()
s_test = test_space.GetVSize()

offsets = mfem.intArray([0, s0, s0+s1])
offsets_test = mfem.intArray([0, s_test])

print('\n'.join(["nNumber of Unknowns",
                 " Trial space,     X0   : " + str(s0) +
                 " (order " + str(trial_order) + ")",
                 " Interface space, Xhat : " + str(s1) +
                 " (order " + str(trace_order) + ")",
                 " Test space,      Y    : " + str(s_test) +
                 " (order " + str(test_order) + ")"]))
x = mfem.BlockVector(offsets)
b = mfem.BlockVector(offsets)
x.Assign(0.0)

# 6. Set up the linear form F(.) which corresponds to the right-hand side of
#    the FEM linear system, which in this case is (f,phi_i) where f=1.0 and
#    phi_i are the basis functions in the test finite element fespace.
one = mfem.ConstantCoefficient(1.0)
F = mfem.LinearForm(test_space)
F.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
F.Assemble()

# 7. Set up the mixed bilinear form for the primal trial unknowns, B0,
#    the mixed bilinear form for the interfacial unknowns, Bhat,
#    the inverse stiffness matrix on the discontinuous test space, Sinv,
#    and the stiffness matrix on the continuous trial space, S0.
ess_bdr = mfem.intArray(mesh.bdr_attributes.Max())
ess_bdr.Assign(1)

B0 = mfem.MixedBilinearForm(x0_space, test_space)
B0.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
B0.Assemble()
B0.EliminateTrialDofs(ess_bdr, x.GetBlock(x0_var), F)
B0.Finalize()

Bhat = mfem.MixedBilinearForm(xhat_space, test_space)
Bhat.AddTraceFaceIntegrator(mfem.TraceJumpIntegrator())
Bhat.Assemble()
Bhat.Finalize()

Sinv = mfem.BilinearForm(test_space)
Sum = mfem.SumIntegrator()
Sum.AddIntegrator(mfem.DiffusionIntegrator(one))
Sum.AddIntegrator(mfem.MassIntegrator(one))

Sinv.AddDomainIntegrator(mfem.InverseIntegrator(Sum))
Sinv.Assemble()
Sinv.Finalize()

S0 = mfem.BilinearForm(x0_space)
S0.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
S0.Assemble()
S0.EliminateEssentialBC(ess_bdr)
S0.Finalize()

matB0 = B0.SpMat()
matBhat = Bhat.SpMat()
matSinv = Sinv.SpMat()
matS0 = S0.SpMat()

# 8. Set up the 1x2 block Least Squares DPG operator, B = [B0  Bhat],
#    the normal equation operator, A = B^t Sinv B, and
#    the normal equation right-hand-size, b = B^t Sinv F.
B = mfem.BlockOperator(offsets_test, offsets)
B.SetBlock(0, 0, matB0)
B.SetBlock(0, 1, matBhat)
A = mfem.RAPOperator(B, matSinv, B)

SinvF = mfem.Vector(s_test)
matSinv.Mult(F, SinvF)
B.MultTranspose(SinvF, b)


# 9. Set up a block-diagonal preconditioner for the 2x2 normal equation
#
#        [ S0^{-1}     0     ]
#        [   0     Shat^{-1} ]      Shat = (Bhat^T Sinv Bhat)
#
#    corresponding to the primal (x0) and interfacial (xhat) unknowns.

# RAP_R replaces RAP, since RAP has two definition one accept
# pointer and the other accept reference. From Python, two
# can not be distingished..
Shat = mfem.RAP_R(matBhat, matSinv, matBhat)
prec_rtol = 1e-3
prec_maxit = 200
S0inv = mfem.CGSolver()
S0inv.SetOperator(matS0)
S0inv.SetPrintLevel(-1)
S0inv.SetRelTol(prec_rtol)
S0inv.SetMaxIter(prec_maxit)

Shatinv = mfem.CGSolver()
Shatinv.SetOperator(Shat)
Shatinv.SetPrintLevel(-1)
Shatinv.SetRelTol(prec_rtol)
Shatinv.SetMaxIter(prec_maxit)
# Disable 'iterative_mode' when using CGSolver (or any IterativeSolver) as
# a preconditioner:
S0inv.iterative_mode = False
Shatinv.iterative_mode = False

P = mfem.BlockDiagonalPreconditioner(offsets)
P.SetDiagonalBlock(0, S0inv)
P.SetDiagonalBlock(1, Shatinv)

# 10. Solve the normal equation system using the PCG iterative solver.
#     Check the weighted norm of residual for the DPG least square problem.
#     Wrap the primal variable in a GridFunction for visualization purposes.
mfem.PCG(A, P, b, x, 1, 200, 1e-12, 0.0)

LSres = mfem.Vector(s_test)
B.Mult(x, LSres)
LSres -= F
res = sqrt(matSinv.InnerProduct(LSres, LSres))

print("|| B0*x0 + Bhat*xhat - F ||_{S^-1} = " + str(res))

x0 = mfem.GridFunction()
x0.MakeRef(x0_space, x.GetBlock(x0_var), 0)

#  11. Save the refined mesh and the solution. This output can be viewed
#     later using GLVis: "glvis -m refined.mesh -g sol.gf".
mesh.Print('refined.mesh', 8)
x0.Save('sol.gf', 8)

if (visualization):
    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock.precision(8)
    sol_sock.send_solution(mesh,  x0)
