'''
   MFEM example 6
      This is a version of Example 1 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np

order = 1
path = dirname((__file__))
meshfile = expanduser(join(path, '..', 'data', 'star.mesh'))

mesh = mfem.Mesh(meshfile, 1, 1)

dim = mesh.Dimension()
sdim = mesh.SpaceDimension()

# 3. Since a NURBS mesh can currently only be refined uniformly, we need to
#    convert it to a piecewise-polynomial curved mesh. First we refine the
#    NURBS mesh a bit more and then project the curvature to quadratic Nodes.
if (mesh.NURBSext):
    for i in range(2):
        mesh.UniformRefinement()
    mesh.SetCurvature(2)

# 4. Define a finite element space on the mesh. The polynomial order is
#      one (linear) by default, but this can be changed on the command line.
fec = mfem.H1_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fec)

# 5. As in Example 1, we set up bilinear and linear forms corresponding to
#    the Laplace problem -\Delta u = 1. We don't assemble the discrete
#    problem yet, this will be done in the main loop.
a = mfem.BilinearForm(fespace)
b = mfem.LinearForm(fespace)

one = mfem.ConstantCoefficient(1.0)
zero = mfem.ConstantCoefficient(0.0)

integ = mfem.DiffusionIntegrator(one)
a.AddDomainIntegrator(integ)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))

# 6. The solution vector x and the associated finite element grid function
#    will be maintained over the AMR iterations. We initialize it to zero.
x = mfem.GridFunction(fespace)
x.Assign(0.0)

#  7. All boundary attributes will be used for essential (Dirichlet) BC.
ess_bdr = intArray(mesh.bdr_attributes.Max())
ess_bdr.Assign(1)


# 9. Set up an error estimator. Here we use the Zienkiewicz-Zhu estimator
#    that uses the ComputeElementFlux method of the DiffusionIntegrator to
#    recover a smoothed flux (gradient) that is subtracted from the element
#    flux to get an error indicator. We need to supply the space for the
#    smoothed flux: an (H1)^sdim (i.e., vector-valued) space is used here.
flux_fespace = mfem.FiniteElementSpace(mesh, fec, sdim)

# own_flux_fes = False indicate flux_fespace is passed by reference
# this is actually default action, but for the sake of explanaiton
# it is explicitly set. If you want to pass pointer use own_flux_fes = True
estimator = mfem.ZienkiewiczZhuEstimator(integ, x, flux_fespace,
                                         own_flux_fes=False)
estimator.SetAnisotropic()

# 10. A refiner selects and refines elements based on a refinement strategy.
#     The strategy here is to refine elements with errors larger than a
#     fraction of the maximum element error. Other strategies are possible.
#     The refiner will call the given error estimator.
refiner = mfem.ThresholdRefiner(estimator)
refiner.SetTotalErrorFraction(0.7)


# 11. The main AMR loop. In each iteration we solve the problem on the
#     current mesh, visualize the solution, and refine the mesh.
max_dofs = 50000
it = 0
while True:
    cdofs = fespace.GetTrueVSize()
    print("AMR iteration " + str(it))
    print("Number of unknowns: " + str(cdofs))

    # 12. Assemble the stiffness matrix and the right-hand side.
    a.Assemble()
    b.Assemble()

    # 13. Set Dirichlet boundary values in the GridFunction x.
    #     Determine the list of Dirichlet true DOFs in the linear system.
    ess_tdof_list = intArray()
    x.ProjectBdrCoefficient(zero, ess_bdr)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 14. Create the linear system: eliminate boundary conditions, constrain
    #     hanging nodes and possibly apply other transformations. The system
    #     will be solved for true (unconstrained) DOFs only.
    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    copy_interior = 1
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior)

    #  15. Define a simple symmetric Gauss-Seidel preconditioner and use it to
    #     solve the linear system with PCG.
    AA = mfem.OperatorHandle2SparseMatrix(A)
    M = mfem.GSSmoother(AA)
    mfem.PCG(AA, M, B, X, 3, 200, 1e-12, 0.0)

    # 16. After solving the linear system, reconstruct the solution as a
    #     finite element GridFunction. Constrained nodes are interpolated
    #     from true DOFs (it may therefore happen that x.Size() >= X.Size()).
    a.RecoverFEMSolution(X, b, x)
    if (cdofs > max_dofs):
        print("Reached the maximum number of dofs. Stop.")
        break
    # 18. Call the refiner to modify the mesh. The refiner calls the error
    #     estimator to obtain element errors, then it selects elements to be
    #     refined and finally it modifies the mesh. The Stop() method can be
    #     used to determine if a stopping criterion was met.
    refiner.Apply(mesh)
    if (refiner.Stop()):
        print("Stopping criterion satisfied. Stop.")
        break
    # 19. Update the space to reflect the new state of the mesh. Also,
    #     interpolate the solution x so that it lies in the new space but
    #     represents the same function. This saves solver iterations later
    #     since we'll have a good initial guess of x in the next step.
    #     Internally, FiniteElementSpace::Update() calculates an
    #     interpolation matrix which is then used by GridFunction::Update().
    fespace.Update()
    x.Update()

    # 20. Inform also the bilinear and linear forms that the space has
    #     changed.
    a.Update()
    b.Update()

    it = it + 1
