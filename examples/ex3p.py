'''
   MFEM example 3 parallel version

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
from mpi4py import MPI
import numpy as np
from numpy import sin, array

freq = 1.0
kappa = np.pi*freq

static_cond = False
order = 1

num_proc = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)
verbose = (myid == 0)


#   3. Read the (serial) mesh from the given mesh file on all processors.  We
#      can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
#      and volume meshes with the same code.
meshfile = expanduser(join(dirname(__file__), '..', 'data', 'beam-tet.mesh'))
mesh = mfem.Mesh(meshfile, 1, 1)

dim = mesh.Dimension()
sdim = mesh.SpaceDimension()


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

#  4. Refine the serial mesh on all processors to increase the resolution. In
#     this example we do 'ref_levels' of uniform refinement. We choose
#     'ref_levels' to be the largest number that gives a final mesh with no
#     more than 1,000 elements.


ref_levels = int(np.floor(np.log(1000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

#   5. Define a parallel mesh by a partitioning of the serial mesh. Refine
#      this mesh further in parallel to increase the resolution. Once the
#      parallel mesh is defined, the serial mesh can be deleted. Tetrahedral
#      meshes need to be reoriented before we can define high-order Nedelec
#      spaces on them.

pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
par_ref_levels = 2
for l in range(par_ref_levels):
    pmesh.UniformRefinement()
pmesh.ReorientTetMesh()

#  6. Define a finite element space on the mesh. Here we use the Nedelec
#     finite elements of the specified order.

fec = mfem.ND_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)
size = fespace.GlobalTrueVSize()

if verbose:  # note that size should be evaulated on all nodes
    print("Number of finite element unknowns: " + str(size))

# 7. Determine the list of true (i.e. conforming) essential boundary dofs.
#    In this example, the boundary conditions are defined by marking all
#    the boundary attributes from the mesh as essential (Dirichlet) and
#    converting them to a list of true dofs.

ess_tdof_list = intArray()
if mesh.bdr_attributes.Size():
    ess_bdr = intArray(mesh.bdr_attributes.Max())
    ess_bdr.Assign(1)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

# 8. Set up the linear form b(.) which corresponds to the right-hand side
#    of the FEM linear system, which in this case is (f,phi_i) where f is
#    given by the function f_exact and phi_i are the basis functions in the
#    finite element fespace.

b = mfem.ParLinearForm(fespace)
f = f_exact()
dd = mfem.VectorFEDomainLFIntegrator(f)
b.AddDomainIntegrator(dd)
b.Assemble()

# 9. Define the solution vector x as a finite element grid function
#    corresponding to fespace. Initialize x by projecting the exact
#    solution. Note that only values from the boundary edges will be used
#    when eliminating the non-homogeneous boundary condition to modify the
#    r.h.s. vector b.

x = mfem.ParGridFunction(fespace)
E = E_exact()
x.ProjectCoefficient(E)

# 10. Set up the bilinear form corresponding to the EM diffusion operator
#       curl muinv curl + sigma I, by adding the curl-curl and the mass domain
#       integrators.

muinv = mfem.ConstantCoefficient(1.0)
sigma = mfem.ConstantCoefficient(1.0)
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))

# 11. Assemble the bilinear form and the corresponding linear system,
#       applying any necessary transformations such as: eliminating boundary
#       conditions, applying conforming constraints for non-conforming AMR,
#       static condensation, etc.
if (static_cond):
    a.EnableStaticCondensation()
a.Assemble()

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

if verbose:
    print("Size of linear system: " + str(A.GetGlobalNumRows()))

# 12. Define and apply a parallel PCG solver for AX=B with the AMS
#     preconditioner from hypre.
prec_fespace = (a.SCParFESpace() if a.StaticCondensationIsEnabled()
                else fespace)
ams = mfem.HypreAMS(A, prec_fespace)
pcg = mfem.HyprePCG(A)
pcg.SetTol(1e-12)
pcg.SetMaxIter(500)
pcg.SetPrintLevel(2)
pcg.SetPreconditioner(ams)
pcg.Mult(B, X)

# 13. Recover the parallel grid function corresponding to X. This is the
#     local finite element solution on each processor.

a.RecoverFEMSolution(X, b, x)

# 12. Compute and print the L^2 norm of the error.
err = x.ComputeL2Error(E)
if verbose:  # note that err should be evaulated on all nodes
    print("|| E_h - E ||_{L^2} = " + "{:g}".format(err))

x.Save('sol.'+smyid)
pmesh.Print('mesh.'+smyid)
