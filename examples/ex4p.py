'''
   MFEM example 4p

   See c++ version in the MFEM library for more detail 

   How to run:
      mpirun -np 2 python <arguments>

   Example of arguments:
      ex4p.py -m square-disc.mesh
      ex4p.py -m star.mesh
      ex4p.py -m beam-tet.mesh
      ex4p.py -m beam-hex.mesh
      ex4p.py -m escher.mesh -o 2 -sc
      ex4p.py -m fichera.mesh -o 2 -hb
      ex4p.py -m fichera-q2.vtk
      ex4p.py -m fichera-q3.mesh -o 2 -sc
      ex4p.py -m square-disc-nurbs.mesh -o 3
      ex4p.py -m beam-hex-nurbs.mesh -o 3
      ex4p.py -m periodic-square.mesh -no-bc
      ex4p.py -m periodic-cube.mesh -no-bc
      ex4p.py -m amr-quad.mesh
      ex4p.py -m amr-hex.mesh -o 2 -sc
      ex4p.py -m amr-hex.mesh -o 2 -hb
      ex4p.py -m star-surf.mesh -o 3 -hb
'''
from numpy import sin, array, cos
import numpy as np
from os.path import expanduser, join, dirname
import sys
from mfem.common.arg_parser import ArgParser
import mfem.par as mfem
from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


parser = ArgParser(description='Ex4 ')
parser.add_argument('-m', '--mesh',
                    default='star.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-no-bc', '--dont-impose-bc',
                    action='store_false',  default=True,
                    help="Impose or not essential boundary conditions.")
parser.add_argument("-f", "--frequency",
                    action='store', default=1.0, type=float,
                    help="Set the frequency for the exact solution.")
parser.add_argument("-sc", "--static-condensation",
                    action='store_true', default=False,
                    help="Enable static condensation.")
parser.add_argument("-hb", "--hybridization",
                    action='store_true', default=False,
                    help="Enable hybridization.")
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')
args = parser.parse_args()

order = args.order
set_bc = args.dont_impose_bc
static_cond = args.static_condensation
hybridization = args.hybridization
freq = args.frequency
kappa = np.pi*freq
visualization = args.visualization
if (myid == 0):
    parser.print_options(args)

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

# 3. Read the (serial) mesh from the given mesh file on all processors.  We
#    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
#    and volume, as well as periodic meshes with the same code.
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
sdim = mesh.SpaceDimension()


#  Coefficients
class E_exact(mfem.VectorPyCoefficient):
    def __init__(self, sdim):
        mfem.VectorPyCoefficient.__init__(self, sdim)

    def EvalValue(self, p):
        dim = p.shape[0]
        x = p[0]
        y = p[1]
        F0 = cos(kappa*x)*sin(kappa*y)
        F1 = cos(kappa*y)*sin(kappa*x)
        if dim == 3:
            return (F0, F1, 0.0)
        else:
            return (F0, F1)


class f_exact(mfem.VectorPyCoefficient):
    def __init__(self, sdim):
        mfem.VectorPyCoefficient.__init__(self, sdim)

    def EvalValue(self, p):
        dim = p.shape[0]
        x = p[0]
        y = p[1]
        temp = 1 + 2*kappa*kappa

        F0 = temp * cos(kappa*x)*sin(kappa*y)
        F1 = temp * cos(kappa*y)*sin(kappa*x)
        if dim == 3:
            return (F0, F1, 0.0)
        else:
            return (F0, F1)


# 4. Refine the serial mesh on all processors to increase the resolution. In
#    this example we do 'ref_levels' of uniform refinement. We choose
#    'ref_levels' to be the largest number that gives a final mesh with no
#    more than 1,000 elements.
ref_levels = int(np.floor(np.log(1000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

# 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
#    this mesh further in parallel to increase the resolution. Once the
#    parallel mesh is defined, the serial mesh can be deleted. Tetrahedral
#    meshes need to be reoriented before we can define high-order Nedelec
#    spaces on them (this is needed in the ADS solver below).
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
par_ref_levels = 2
for l in range(par_ref_levels):
    pmesh.UniformRefinement()
pmesh.ReorientTetMesh()

# 6. Define a parallel finite element space on the parallel mesh. Here we
#    use the Raviart-Thomas finite elements of the specified order.
fec = mfem.RT_FECollection(order-1, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)
size = fespace.GlobalTrueVSize()
if myid == 0:
    print("Number of finite element unknows : " + str(size))

# 7. Determine the list of true (i.e. parallel conforming) essential
#    boundary dofs. In this example, the boundary conditions are defined
#    by marking all the boundary attributes from the mesh as essential
#    (Dirichlet) and converting them to a list of true dofs.
ess_tdof_list = mfem.intArray()
if pmesh.bdr_attributes.Size():
    ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
    if set_bc:
        ess_bdr.Assign(1)
    else:
        ess_bdr.Assign(0)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

# 8. Set up the parallel linear form b(.) which corresponds to the
#    right-hand side of the FEM linear system, which in this case is
#    (f,phi_i) where f is given by the function f_exact and phi_i are the
#    basis functions in the finite element fespace.
f = f_exact(sdim)
b = mfem.ParLinearForm(fespace)
b.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(f))
b.Assemble()

# 9. Define the solution vector x as a parallel finite element grid function
#    corresponding to fespace. Initialize x by projecting the exact
#    solution. Note that only values from the boundary faces will be used
#    when eliminating the non-homogeneous boundary condition to modify the
#    r.h.s. vector b.
x = mfem.ParGridFunction(fespace)
F = E_exact(sdim)
x.ProjectCoefficient(F)

# 10. Set up the parallel bilinear form corresponding to the H(div)
#     diffusion operator grad alpha div + beta I, by adding the div-div and
#     the mass domain integrators.

alpha = mfem.ConstantCoefficient(1.0)
beta = mfem.ConstantCoefficient(1.0)
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DivDivIntegrator(alpha))
a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(beta))

# 11. Assemble the parallel bilinear form and the corresponding linear
#     system, applying any necessary transformations such as: parallel
#     assembly, eliminating boundary conditions, applying conforming
#     constraints for non-conforming AMR, static condensation,
#     hybridization, etc.
if (static_cond):
    a.EnableStaticCondensation()
elif (hybridization):
    hfec = mfem.DG_Interface_FECollection(order-1, dim)
    hfes = mfem.ParFiniteElementSpace(mesh, hfec)
    a.EnableHybridization(hfes, mfem.NormalTraceJumpIntegrator(),
                          ess_tdof_list)
a.Assemble()

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

glob_size = A.GetGlobalNumRows()
if myid == 0:
    print("Size of linear system: " + str(glob_size))

# 12. Define and apply a parallel PCG solver for A X = B with the 2D AMS or
#     the 3D ADS preconditioners from hypre. If using hybridization, the
#     system is preconditioned with hypre's BoomerAMG.

pcg = mfem.CGSolver(A.GetComm())
pcg.SetOperator(A)
pcg.SetRelTol(1e-12)
pcg.SetMaxIter(500)
pcg.SetPrintLevel(1)
if hybridization:
    prec = mfem.HypreBoomerAMG(A)
else:
    if a.StaticCondensationIsEnabled():
        prec_fespace = a.SCParFESpace()
    else:
        prec_fespace = fespace
    if dim == 2:
        prec = mfem.HypreAMS(A, prec_fespace)
    else:
        prec = mfem.HypreADS(A, prec_fespace)

pcg.SetPreconditioner(prec)
pcg.Mult(B, X)

# 13. Recover the parallel grid function corresponding to X. This is the
#     local finite element solution on each processor.
a.RecoverFEMSolution(X, b, x)

# 14. Compute and print the L^2 norm of the error.
err = x.ComputeL2Error(F)
if myid == 0:
    print("|| F_h - F ||_{L^2} = " + str(err))

# 15. Save the refined mesh and the solution in parallel. This output can
#     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
sol_name = "sol."+smyid
pmesh.Print(mesh_name, 8)
x.Save(sol_name, 8)

# 16. Send the solution by socket to a GLVis server.
if visualization:
    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock.send_text("parallel " + str(num_procs) + " " + str(myid))
    sol_sock.precision(8)
    sol_sock.send_solution(pmesh,   x)
