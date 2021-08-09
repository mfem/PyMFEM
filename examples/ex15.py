'''
   MFEM example 15

   How to run:
      python <arguments>

   Example of arguments:
      ex15.py
      ex15.py -o 1 -y 0.4
      ex15.py -o 4 -y 0.1
      ex15.py -n 5
      ex15.py -p 1 -n 3

      Other meshes:

      ex15.py -m square-disc-nurbs.mesh
      ex15.py -m disc-nurbs.mesh
      ex15.py -m fichera.mesh -tf 0.3
      ex15.py -m ball-nurbs.mesh -tf 0.3
      ex15.py -m mobius-strip.mesh
      ex15.py -m amr-quad.mesh

      Conforming meshes (no derefinement):

      ex15.py -m square-disc.mesh
      ex15.py -m escher.mesh -r 2 -tf 0.3

'''

import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from numpy import cos, sin, pi, exp, sqrt, arctan

import mfem.ser as mfem
from mfem.ser import intArray

parser = ArgParser(description='Ex15')
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
                    action='store', default=5e-3, type=float,
                    help="Maximum element error")
parser.add_argument("-y", "--hysteresis",
                    action='store', default=0.15, type=float,
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
parser.print_options(args)


alpha = 0.02


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


def UpdateProblem(mesh, fespace, x, a, b):
    # Update the space: recalculate the number of DOFs and construct a matrix
    # that will adjust any GridFunctions to the new mesh state.
    fespace.Update()

    # Interpolate the solution on the new mesh by applying the transformation
    # matrix computed in the finite element space. Multiple GridFunctions could
    # be updated here.
    x.Update()

    # Free any transformation matrices to save memory.
    fespace.UpdatesFinished()

    # Inform the linear and bilinear forms that the space has changed.
    a.Update()
    b.Update()


# 2. Read the mesh from the given mesh file on all processors. We can handle
#    triangular, quadrilateral, tetrahedral, hexahedral, surface and volume
#    meshes with the same code
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
sdim = mesh.SpaceDimension()

# 3. Project a NURBS mesh to a piecewise-quadratic curved mesh. Make sure
#    that the mesh is non-conforming if it has quads or hexes and refine it
if (mesh.NURBSext):
    mesh.UniformRefinement()
    if ref_levels > 0:
        ref_levels = ref_levels-1
    mesh.SetCurvature(2)

mesh.EnsureNCMesh()
for l in range(ref_levels):
    mesh.UniformRefinement()

# 4. All boundary attributes will be used for essential (Dirichlet) BC
ess_bdr = intArray(mesh.bdr_attributes.Max())
ess_bdr.Assign(1)

# 5. Define a finite element space on the mesh. The polynomial order is one
#   (linear) by default, but this can be changed on the command line.
fec = mfem.H1_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fec)

# 6. As in Example 1p, we set up bilinear and linear forms corresponding to
#    the Laplace problem -\Delta u = 1. We don't assemble the discrete
#    problem yet, this will be done in the inner loop.
a = mfem.BilinearForm(fespace)
b = mfem.LinearForm(fespace)

one = mfem.ConstantCoefficient(1.0)
bdr = BdrCoefficient()
rhs = RhsCoefficient()

integ = mfem.DiffusionIntegrator(one)
a.AddDomainIntegrator(integ)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(rhs))

# 7. The solution vector x and the associated finite element grid function
#    will be maintained over the AMR iterations.
x = mfem.GridFunction(fespace)

# 8. Connect to GLVis.
if visualization:
    sout = mfem.socketstream("localhost", 19916)
    sout.precision(8)

# 9. As in Example 6, we set up a Zienkiewicz-Zhu estimator that will be
#    used to obtain element error indicators. The integrator needs to
#    provide the method ComputeElementFlux. The smoothed flux space is a
#    vector valued H1 space here.
flux_fespace = mfem.FiniteElementSpace(mesh, fec, sdim)

# own_flux_fes = False indicate flux_fespace is passed by reference
# this is actually default action, but for the sake of explanaiton
# it is explicitly set. If you want to pass pointer use own_flux_fes = True
estimator = mfem.ZienkiewiczZhuEstimator(integ, x, flux_fespace,
                                         own_flux_fes=False)

# 10. As in Example 6, we also need a refiner. This time the refinement
#     strategy is based on a fixed threshold that is applied locally to each
#     element. The global threshold is turned off by setting the total error
#     fraction to zero. We also enforce a maximum refinement ratio between
#     adjacent elements.
refiner = mfem.ThresholdRefiner(estimator)
refiner.SetTotalErrorFraction(0.0)
refiner.SetLocalErrorGoal(max_elem_error)
refiner.PreferConformingRefinement()
refiner.SetNCLimit(nc_limit)

# 11. A derefiner selects groups of elements that can be coarsened to form
#     a larger element. A conservative enough threshold needs to be set to
#     prevent derefining elements that would immediately be refined again.
derefiner = mfem.ThresholdDerefiner(estimator)
derefiner.SetThreshold(hysteresis * max_elem_error)
derefiner.SetNCLimit(nc_limit)

# 12. The outer time loop. In each iteration we update the right hand side,
#     solve the problem on the current mesh, visualize the solution and
#     refine the mesh as many times as necessary. Then we derefine any
#     elements which have very small errors.
x.Assign(0.0)
time = 0.0
while (time < t_final + 1e-10):
    print("Time " + str(time) + "\n\nRefinement:")
    bdr.SetTime(time)
    rhs.SetTime(time)

    # Make sure errors will be recomputed in the following.
    refiner.Reset()
    derefiner.Reset()

    # 13. The inner refinement loop. At the end we want to have the current
    #     time step resolved to the prescribed tolerance in each element.
    ref_it = 0
    while(True):
        ref_it = ref_it + 1
        print("Iteration: " + str(ref_it) + ", number of unknowns: "
              + str(fespace.GetVSize()))

        # 14. Recompute the field on the current mesh: assemble the stiffness
        #     matrix and the right-hand side.
        a.Assemble()
        b.Assemble()

        # 15. Project the exact solution to the essential boundary DOFs.
        x.ProjectBdrCoefficient(bdr, ess_bdr)
        # 16. Create and solve the linear system.
        ess_tdof_list = intArray()
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

        A = mfem.OperatorPtr()
        B = mfem.Vector()
        X = mfem.Vector()
        a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

        AA = mfem.OperatorHandle2SparseMatrix(A)
        M = mfem.GSSmoother(AA)
        mfem.PCG(AA, M, B, X, 0, 500, 1e-12, 0.0)

        # 17. Extract the local solution on each processor.
        a.RecoverFEMSolution(X, b, x)

        # 18. Send the solution by socket to a GLVis server
        if visualization:
            sout.precision(8)
            sout.send_solution(mesh, x)

        # 19. Apply the refiner on the mesh. The refiner calls the error
        #     estimator to obtain element errors, then it selects elements to
        #     be refined and finally it modifies the mesh. The Stop() method
        #     determines if all elements satisfy the local threshold.
        refiner.Apply(mesh)
        if refiner.Stop():
            break

        # 20. Update the space and interpolate the solution.
        UpdateProblem(mesh, fespace, x, a, b)
    # 21. Use error estimates from the last inner iteration to check for
    #     possible derefinements. The derefiner works similarly as the
    #     refiner. The errors are not recomputed because the mesh did not
    #     change (and also the estimator was not Reset() at this time).
    if derefiner.Apply(mesh):
        print("Derefined elements.")
        #  22. Update the space and interpolate the solution.
        UpdateProblem(mesh, fespace, x, a, b)
    a.Update()
    b.Update()
    time = time + 0.01
