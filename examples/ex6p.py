'''
   MFEM example 6
      This is a version of Example 1 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
from mpi4py import MPI
import numpy as np

num_proc = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)
verbose = (myid == 0)

order = 1
meshfile = expanduser(join(dirname(__file__), '..', 'data', 'star.mesh'))
mesh = mfem.Mesh(meshfile, 1, 1)

dim = mesh.Dimension()
sdim = mesh.SpaceDimension()

# 4. Refine the serial mesh on all processors to increase the resolution.
#    Also project a NURBS mesh to a piecewise-quadratic curved mesh. Make
#    sure that the mesh is non-conforming.

if (mesh.NURBSext):
    mesh.UniformRefinement()
    mesh.SetCurvature(2)
mesh.EnsureNCMesh()

# 5. Define a parallel mesh by partitioning the serial mesh.
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
ess_bdr = intArray(pmesh.bdr_attributes.Max())
ess_bdr.Assign(1)

# 6. Define a finite element space on the mesh. The polynomial order is
#    one (linear) by default, but this can be changed on the command line.
fec = mfem.H1_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)

# 7. As in Example 1p, we set up bilinear and linear forms corresponding to
#    the Laplace problem -\Delta u = 1. We don't assemble the discrete
#    problem yet, this will be done in the main loop.
a = mfem.ParBilinearForm(fespace)
b = mfem.ParLinearForm(fespace)

one = mfem.ConstantCoefficient(1.0)

integ = mfem.DiffusionIntegrator(one)
a.AddDomainIntegrator(integ)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))

# 8. The solution vector x and the associated finite element grid function
#    will be maintained over the AMR iterations. We initialize it to zero.
x = mfem.ParGridFunction(fespace)
x.Assign(0.0)

# 10. Set up an error estimator. Here we use the Zienkiewicz-Zhu estimator
#     with L2 projection in the smoothing step to better handle hanging
#     nodes and parallel partitioning. We need to supply a space for the
#     discontinuous flux (L2) and a space for the smoothed flux (H(div) is
#     used here).
flux_fec = mfem.L2_FECollection(order, dim)
flux_fes = mfem.ParFiniteElementSpace(pmesh, flux_fec, sdim)
smooth_flux_fec = mfem.RT_FECollection(order-1, dim)
smooth_flux_fes = mfem.ParFiniteElementSpace(pmesh, smooth_flux_fec)


estimator = mfem.L2ZienkiewiczZhuEstimator(integ, x, flux_fes,
                                           smooth_flux_fes)
# 11. A refiner selects and refines elements based on a refinement strategy.
#     The strategy here is to refine elements with errors larger than a
#     fraction of the maximum element error. Other strategies are possible.
#     The refiner will call the given error estimator.
refiner = mfem.ThresholdRefiner(estimator)
refiner.SetTotalErrorFraction(0.7)

# 11. The main AMR loop. In each iteration we solve the problem on the
#     current mesh, visualize the solution, and refine the mesh.
max_dofs = 100000
it = 0
while True:
    global_dofs = fespace.GlobalTrueVSize()
    if (myid == 0):
        print("AMR iteration " + str(it))
        print("Number of unknowns: " + str(global_dofs))

    # 12. Assemble the stiffness matrix and the right-hand side.
    a.Assemble()
    b.Assemble()

    # 13. Set Dirichlet boundary values in the GridFunction x.
    #     Determine the list of Dirichlet true DOFs in the linear system.
    ess_tdof_list = intArray()
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 14. Create the linear system: eliminate boundary conditions, constrain
    #     hanging nodes and possibly apply other transformations. The system
    #     will be solved for true (unconstrained) DOFs only.
    A = mfem.HypreParMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    copy_interior = 1
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior)

    #  15. Define and apply a parallel PCG solver for AX=B with the BoomerAMG
    #      preconditioner from hypre.

    amg = mfem.HypreBoomerAMG()
    amg.SetPrintLevel(0)
    pcg = mfem.CGSolver(A.GetComm())
    pcg.SetPreconditioner(amg)
    pcg.SetOperator(A)
    pcg.SetRelTol(1e-6)
    pcg.SetMaxIter(200)
    pcg.SetPrintLevel(3)
    pcg.Mult(B, X)

    # 16. Extract the parallel grid function corresponding to the finite element
    #     approximation X. This is the local solution on each processor.
    a.RecoverFEMSolution(X, b, x)
    if (global_dofs > max_dofs):
        if (myid == 0):
            print("Reached the maximum number of dofs. Stop.")
        break

    refiner.Apply(pmesh)
    if (refiner.Stop()):
        if myid == 0:
            print("Stopping criterion satisfied. Stop.")
        break
    # 19. Update the finite element space (recalculate the number of DOFs,
    #     etc.) and create a grid function update matrix. Apply the matrix
    #     to any GridFunctions over the space. In this case, the update
    #     matrix is an interpolation matrix so the updated GridFunction will
    #     still represent the same function as before refinement.
    fespace.Update()
    x.Update()

    # 20. Load balance the mesh, and update the space and solution. Currently
    #     available only for nonconforming meshes.
    if pmesh.Nonconforming():
        pmesh.Rebalance()
        # Update the space and the GridFunction. This time the update matrix
        # redistributes the GridFunction among the processors.
        fespace.Update()
        x.Update()

    # 21. Inform also the bilinear and linear forms that the space has
    #     changed.
    a.Update()
    b.Update()

    if ((it + 1) % 5) == 0:
        pmesh.ParPrint("ex6p-checkpoint."+smyid, 8)
        if myid == 0:
            print("Checkpoint saved")
    it = it + 1
