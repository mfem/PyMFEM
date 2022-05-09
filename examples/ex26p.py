'''
   MFEM example 26
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(order_refinements=2,
        geometric_refinements=0,
        mesh_file='',
        visualization=True,
        device='cpu'):

    class DiffusionMultigrid(mfem.PyGeometricMultigrid):
        def __init__(self, fespaces, ess_bdr):
            mfem.PyGeometricMultigrid.__init__(self, fespaces)

            self.smoothers = []
            self.ConstructCoarseOperatorAndSolver(fespaces.GetFESpaceAtLevel(0),
                                                  ess_bdr)

            for level in range(1, fespaces.GetNumLevels()):
                self.ConstructOperatorAndSmoother(fespaces.GetFESpaceAtLevel(level),
                                                  ess_bdr)

        def ConstructBilinearForm(self, fespace, ess_bdr, partial_assembly):
            form = mfem.ParBilinearForm(fespace)
            if partial_assembly:
                form.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
            form.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
            form.Assemble()
            self.AppendBilinearForm(form)

            ess = mfem.intArray()
            self.AppendEssentialTDofs(ess)
            fespace.GetEssentialTrueDofs(ess_bdr, ess)

        def ConstructCoarseOperatorAndSolver(self, coarse_fespace,
                                             ess_bdr):

            self.ConstructBilinearForm(coarse_fespace, ess_bdr, False)

            hypreCoarseMat = mfem.HypreParMatrix()
            self.bfs[-1].FormSystemMatrix(self.essentialTrueDofs[-1],
                                          hypreCoarseMat)

            amg = mfem.HypreBoomerAMG(hypreCoarseMat)
            amg.SetPrintLevel(-1)

            pcg = mfem.CGSolver(MPI.COMM_WORLD)
            pcg.SetPrintLevel(-1)
            pcg.SetMaxIter(200)
            pcg.SetRelTol(sqrt(1e-4))
            pcg.SetAbsTol(0.0)
            pcg.SetOperator(hypreCoarseMat)
            pcg.SetPreconditioner(amg)

            self.smoothers.append((pcg, hypreCoarseMat))
            self.AddLevel(hypreCoarseMat, pcg, False, False)
            self.amg = amg

        def ConstructOperatorAndSmoother(self, fespace, ess_bdr):
            self.ConstructBilinearForm(fespace, ess_bdr, True)

            opr = mfem.OperatorPtr()
            opr.SetType(mfem.Operator.ANY_TYPE)

            self.bfs[-1].FormSystemMatrix(self.essentialTrueDofs[-1], opr)
            opr.SetOperatorOwner(False)

            diag = mfem.Vector(fespace.GetTrueVSize())
            self.bfs[-1].AssembleDiagonal(diag)

            smoother = mfem.OperatorChebyshevSmoother(opr.Ptr(), diag,
                                                      self.essentialTrueDofs[-1],
                                                      2,
                                                      fespace.GetParMesh().GetComm())

            self.smoothers.append((smoother, opr))
            self.AddLevel(opr.Ptr(), smoother, False, False)

    device = mfem.Device(device)
    if myid == 0:
        device.Print()

    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()

    ref_levels = int(np.floor(np.log(1000. / mesh.GetNE()) / np.log(2.) / dim))
    for i in range(ref_levels):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    pmesh.UniformRefinement()
    pmesh.UniformRefinement()

    # 5. Define a finite element space hierarchy on the mesh. Here we use
    #    continuous Lagrange finite elements. We start with order 1 on the
    #    coarse level and geometrically refine the spaces by the specified
    #    amount. Afterwards, we increase the order of the finite elements
    #    by a factor of 2 for each additional level.
    fec = mfem.H1_FECollection(1, dim)
    coarse_fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    # note ownM and ownFES are False since they are managed by Python GC.
    fespaces = mfem.ParFiniteElementSpaceHierarchy(
        pmesh, coarse_fespace, False, False)

    collections = []
    collections.append(fec)
    for level in range(geometric_refinements):
        fespaces.AddUniformlyRefinedLevel()
    for level in range(order_refinements):
        collections.append(mfem.H1_FECollection(2**(level+1), dim))
        fespaces.AddOrderRefinedLevel(collections[-1])

    global_dof = fespaces.GetFinestFESpace().GlobalTrueVSize()
    if myid == 0:
        print("Number of finite element unknowns: " + str(global_dof))

    # 6. Set up the linear form b(.) which corresponds to the right-hand side of
    #    the FEM linear system, which in this case is (1,phi_i) where phi_i are
    #    the basis functions in the finite element fespace.
    b = mfem.ParLinearForm(fespaces.GetFinestFESpace())
    one = mfem.ConstantCoefficient(1.0)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
    b.Assemble()

    # 7. Define the solution vector x as a finite element grid function
    #    corresponding to fespace. Initialize x with initial guess of zero,
    #    which satisfies the boundary conditions.
    x = mfem.ParGridFunction(fespaces.GetFinestFESpace())
    x.Assign(0.0)

    # 8. Create the multigrid operator using the previously created
    #    FiniteElementSpaceHierarchy and additional boundary information. This
    #    operator is then used to create the MultigridSolver as a preconditioner
    #    in the iterative solver.
    ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)

    M = DiffusionMultigrid(fespaces, ess_bdr)
    M.SetCycleType(mfem.Multigrid.CycleType_VCYCLE, 1, 1)

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()

    M.FormFineLinearSystem(x, b, A, X, B)
    if myid == 0:
        print("Size of linear system: " + str(A.Height()))

    # 9. Solve the linear system A X = B.
    cg = mfem.CGSolver(MPI.COMM_WORLD)
    cg.SetRelTol(1e-12)
    cg.SetMaxIter(200)
    cg.SetPrintLevel(1)
    cg.SetOperator(A.Ptr())
    cg.SetPreconditioner(M)

    import sys
    sys.setrecursionlimit(200000)
    cg.Mult(B, X)

    # 10. Recover the solution as a finite element grid function.
    M.RecoverFineFEMSolution(X, b, x)

    # 11. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g sol.gf".
    fespaces.GetFinestFESpace().GetMesh().Print("mesh"+smyid, 8)
    x.Save("sol"+smyid, 8)

    # 12. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << fespaces.GetFinestFESpace().GetMesh() << x
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex26p (Multigrid Preconditioner)')

    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-gr', '--geometric-refinement',
                        action='store', default=0, type=int,
                        help="Number of geometric refinements done prior to order refinements.")
    parser.add_argument('-or', '--order-refinement',
                        action='store', default=2, type=int,
                        help="Number of order refinements. Finest level in the hierarchy has order 2^{or}.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    mesh_file = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    run(order_refinements=args.order_refinement,
        geometric_refinements=args.geometric_refinement,
        mesh_file=mesh_file,
        visualization=args.visualization,
        device=args.device)
