'''
   MFEM example 15pe

   How to run:
      mpirun -np 2 python <arguments>

   Example of arguments:
      ex15p.py
      ex15p.py -o 1 -y 0.4
      ex15p.py -o 4 -y 0.1
      ex15p.py -n 5
      ex15p.py -p 1 -n 3

      Other meshes:

      ex15p.py -m square-disc-nurbs.mesh
      ex15p.py -m disc-nurbs.mesh
      ex15p.py -m fichera.mesh -tf 0.5
      ex15p.py -m ball-nurbs.mesh -tf 0.5
      ex15p.py -m mobius-strip.mesh
      ex15p.py -m amr-quad.mesh

      Conforming meshes (no derefinement):

      ex15p.py -m square-disc.mesh
      ex15p.py -m escher.mesh -r 2 -tf 0.3

'''

import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from numpy import cos, sin, pi, exp, sqrt, arctan

import mfem.par as mfem
from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

parser = ArgParser(description='Ex15p')
parser.add_argument('-m', '--mesh',
                    default='star-hilbert.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument("-p", "--problem",
                    action='store', default=0, type=int,
                    help="Problem setup to use: 0 = spherical front, 1 = ball.")
parser.add_argument("-n", "--nfeatures",
                    action='store', default=1, type=int,
                    help="Number of solution features (fronts/balls).")
parser.add_argument('-o', '--order',
                    action='store', default=2, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument("-e", "--max-err",
                    action='store', default=1e-4, type=float,
                    help="Maximum element error")
parser.add_argument("-y", "--hysteresis",
                    action='store', default=0.25, type=float,
                    help="Derefinement safety coefficient.")
parser.add_argument('-r', '--ref-levels',
                    action='store', default=0, type=int,
                    help="Number of inital uniform refinement")
parser.add_argument("-l", "--nc-limit",
                    action='store', default=3, type=int,
                    help="Maximum level of hanging nodes.")
parser.add_argument('-tf', '--t-final',
                    action='store', default=1.0, type=float,
                    help="Final time; start time is 0.")
parser.add_argument('-vis', '--visualization',
                    action='store_true', default=True,
                    help='Enable GLVis visualization')

args = parser.parse_args()
ref_levels = args.ref_levels
problem = args.problem
nfeatures = args.nfeatures
order = args.order
max_elem_error = args.max_err
hysteresis = args.hysteresis
nc_limit = args.nc_limit
t_final = args.t_final
visualization = args.visualization
if myid == 0:
    parser.print_options(args)
    
alpha = 0.02

device = mfem.Device('cpu')
if myid == 0:
    device.Print()

def front(x, y, z, t, dim):
    r = sqrt(x**2 + y**2 + z**2)
    return exp(-0.5 * ((r - t)/alpha) ** 2)


def ball(x, y, z, t, dim):
    r = sqrt(x**2 + y**2 + z**2)
    return -atan(2.*(r - t)/alpha)


def front_laplace(x, y, z, t, dim):
    x2 = x**2
    y2 = y**2
    z2 = z**2
    t2 = t**2
    a2 = alpha**2
    a4 = a2**2
    r = sqrt(x2 + y2 + z2)
    ret = (-exp(- 0.5 * ((r - t)/alpha)**2) / a4 *
           (-2.*t*(x2 + y2 + z2 - (dim-1)*a2/2.)/r + x2 + y2 + z2 + t2 - dim*a2))
    return ret


def ball_laplace(x, y, z, t, dim):
    x2 = x**2
    y2 = y**2
    z2 = z**2
    t2 = t**2
    a2 = alpha**2
    a4 = a2**2
    r = sqrt(x2 + y2 + z2)
    den = (-a2 - 4.*(x2 + y2 + z2 - 2*r*t) - t2)**2

    if dim == 2:
        return 2*alpha*(a2 + t2 - 4*x2 - 4*y2)/r/den
    return 4*alpha*(a2 + t2 - 4.*r*t)/r/den


def composite_func(pt, t, f0, f1):
    dim = len(pt)
    x = pt[0]
    y = pt[1]
    z = 0.0
    if dim == 3:
        z = pt[2]
    if (problem == 0):
        if (nfeatures <= 1):
            return f0(x, y, z, t, dim)
        else:
            i = np.arange(nfeatures)
            x0 = 0.5*cos(2 * pi * i / nfeatures)
            y0 = 0.5*sin(2 * pi * i / nfeatures)
            return np.sum(f0(x - x0, y - y0, z, t, dim))
    else:
        i = np.arange(nfeatures)
        x0 = 0.5*cos(2 * pi * i / nfeatures + pi*t)
        y0 = 0.5*sin(2 * pi * i / nfeatures + pi*t)
        return np.sum(f1(x - x0, y - y0, z, 0.25, dim))


class BdrCoefficient(mfem.PyCoefficientT):
    def EvalValue(self, pt, t):
        return composite_func(pt, t, front, ball)


class RhsCoefficient(mfem.PyCoefficientT):
    def EvalValue(self, pt, t):
        # print 'rhs', composite_func(pt, t, front_laplace, ball_laplace)
        return composite_func(pt, t, front_laplace, ball_laplace)


def UpdateAndRebalance(pmesh, fespace, x, a, b):
    # Update the space: recalculate the number of DOFs and construct a matrix
    # that will adjust any GridFunctions to the new mesh state.
    fespace.Update()

    # Interpolate the solution on the new mesh by applying the transformation
    # matrix computed in the finite element space. Multiple GridFunctions could
    # be updated here.
    x.Update()

    if pmesh.Nonconforming():
        # Load balance the mesh.
        pmesh.Rebalance()
        # Update the space again, this time a GridFunction redistribution matrix
        # is created. Apply it to the solution.
        fespace.Update()
        x.Update()

    # Inform the linear and bilinear forms that the space has changed.
    a.Update()
    b.Update()

    # Free any transformation matrices to save memory.
    fespace.UpdatesFinished()

# 3. Read the (serial) mesh from the given mesh file on all processors.  We
#    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
#    and volume meshes with the same code.


meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
sdim = mesh.SpaceDimension()

# 4. Project a NURBS mesh to a piecewise-quadratic curved mesh. Make sure
#    that the mesh is non-conforming if it has quads or hexes and refine it
if (mesh.NURBSext):
    mesh.UniformRefinement()
    if ref_levels > 0:
        ref_levels = ref_levels-1
    mesh.SetCurvature(2)

mesh.EnsureNCMesh()
for l in range(ref_levels):
    mesh.UniformRefinement()

# 5. Define a parallel mesh by partitioning the serial mesh.  Once the
#    parallel mesh is defined, the serial mesh can be deleted.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
ess_bdr.Assign(1)

# 6. Define a finite element space on the mesh. The polynomial order is one
#   (linear) by default, but this can be changed on the command line.
fec = mfem.H1_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)

# 7. As in Example 1p, we set up bilinear and linear forms corresponding to
#    the Laplace problem -\Delta u = 1. We don't assemble the discrete
#    problem yet, this will be done in the inner loop.
a = mfem.ParBilinearForm(fespace)
b = mfem.ParLinearForm(fespace)

one = mfem.ConstantCoefficient(1.0)
bdr = BdrCoefficient()
rhs = RhsCoefficient()

integ = mfem.DiffusionIntegrator(one)
a.AddDomainIntegrator(integ)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(rhs))

# 8. The solution vector x and the associated finite element grid function
#    will be maintained over the AMR iterations.
x = mfem.ParGridFunction(fespace)

# 9. Connect to GLVis.
if visualization:
    sout = mfem.socketstream("localhost", 19916)
    sout.precision(8)

# 10. As in Example 6p, we set up a Zienkiewicz-Zhu estimator that will be
#     used to obtain element error indicators. The integrator needs to
#     provide the method ComputeElementFlux. We supply an L2 space for the
#     discontinuous flux and an H(div) space for the smoothed flux.
flux_fec = mfem.L2_FECollection(order, dim)
flux_fes = mfem.ParFiniteElementSpace(pmesh, flux_fec, sdim)
smooth_flux_fec = mfem.RT_FECollection(order-1, dim)
smooth_flux_fes = mfem.ParFiniteElementSpace(pmesh, smooth_flux_fec)
estimator = mfem.L2ZienkiewiczZhuEstimator(integ, x, flux_fes,
                                           smooth_flux_fes)

# 11. As in Example 6p, we also need a refiner. This time the refinement
#     strategy is based on a fixed threshold that is applied locally to each
#     element. The global threshold is turned off by setting the total error
#     fraction to zero. We also enforce a maximum refinement ratio between
#     adjacent elements.
refiner = mfem.ThresholdRefiner(estimator)
refiner.SetTotalErrorFraction(0.0)   # purely local threshold
refiner.SetLocalErrorGoal(max_elem_error)
refiner.PreferConformingRefinement()
refiner.SetNCLimit(nc_limit)

# 12. A derefiner selects groups of elements that can be coarsened to form
#     a larger element. A conservative enough threshold needs to be set to
#      prevent derefining elements that would immediately be refined again.
derefiner = mfem.ThresholdDerefiner(estimator)
derefiner.SetThreshold(hysteresis * max_elem_error)
derefiner.SetNCLimit(nc_limit)

# 13. The outer time loop. In each iteration we update the right hand side,
#     solve the problem on the current mesh, visualize the solution and
#     refine the mesh as many times as necessary. Then we derefine any
#     elements which have very small errors.
# x.Assign(0.0)
time = 0.0
while (time < t_final + 1e-10):
    if myid == 0:
        print("Time " + str(time) + "\n\nRefinement:")
    bdr.SetTime(time)
    rhs.SetTime(time)

    # Make sure errors will be recomputed in the following.
    refiner.Reset()
    derefiner.Reset()

    # 14. The inner refinement loop. At the end we want to have the current
    #     time step resolved to the prescribed tolerance in each element.
    ref_it = 1
    while(True):
        global_dofs = fespace.GlobalTrueVSize()
        if myid == 0:
            print_txt = ("Iteration: " + str(ref_it) + ", number of unknowns: "
                         + str(global_dofs))

        # 15. Recompute the field on the current mesh: assemble the stiffness
        #     matrix and the right-hand side.
        a.Assemble()
        b.Assemble()

        # 16. Project the exact solution to the essential DOFs.
        x.ProjectBdrCoefficient(bdr, ess_bdr)

        # 17. Create and solve the parallel linear system.
        ess_tdof_list = mfem.intArray()
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

        A = mfem.HypreParMatrix()
        B = mfem.Vector()
        X = mfem.Vector()
        a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

        amg = mfem.HypreBoomerAMG(A)
        amg.SetPrintLevel(0)
        pcg = mfem.HyprePCG(A)
        pcg.SetTol(1e-12)
        pcg.SetMaxIter(200)
        pcg.SetPrintLevel(0)
        pcg.SetPreconditioner(amg)
        pcg.Mult(B, X)

        # 18. Extract the local solution on each processor.
        a.RecoverFEMSolution(X, b, x)

        # 19. Send the solution by socket to a GLVis server
        if visualization:
            sout.send_text("parallel " + str(num_procs) + " " + str(myid))
            sout.precision(8)
            sout.send_solution(pmesh, x)

        # 20. Apply the refiner on the mesh. The refiner calls the error
        #     estimator to obtain element errors, then it selects elements to
        #     be refined and finally it modifies the mesh. The Stop() method
        #     determines if all elements satisfy the local threshold.
        refiner.Apply(pmesh)
        if myid == 0:
            print(print_txt + ", total error: " +
                  str(estimator.GetTotalError()))
        # 21. Quit the AMR loop if the termination criterion has been met
        if refiner.Stop():
            a.Update()
            break

        # 22. Update the space, interpolate the solution, rebalance the mesh.
        UpdateAndRebalance(pmesh, fespace, x, a, b)
        ref_it = ref_it + 1
    # 21. Use error estimates from the last inner iteration to check for
    #     possible derefinements. The derefiner works similarly as the
    #     refiner. The errors are not recomputed because the mesh did not
    #     change (and also the estimator was not Reset() at this time).
    if derefiner.Apply(pmesh):
        if myid == 0:
            print("Derefined elements.")
        #  22. Update the space and interpolate the solution.
        UpdateAndRebalance(pmesh, fespace, x, a, b)

    time = time + 0.01
