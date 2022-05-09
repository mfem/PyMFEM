'''
   MFEM example 5p

   See c++ version in the MFEM library for more detail 
'''
import time
from numpy import sin, cos, exp
import numpy as np
from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from mfem.common.mpi_debug import nicePrint

from os.path import expanduser, join, dirname

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)
verbose = (myid == 0)


# time.clock deprecated and removed in PY3.8
try:
    clock = time.process_time
except AttributeError:
    clock = time.clock


def run(order=1,
        refine=-1,
        meshfile='',
        visualization=False,
        par_format=False,
        device='cpu',
        numba=False,
        pa=False):

    def pFunc_exact(x):
        xi = float(x[0])
        yi = float(x[1])
        zi = 0.0
        if len(x) == 3:
            zi = x[2]
        from numpy import sin, cos, exp
        return exp(xi)*sin(yi)*cos(zi)

    class uFunc_ex(mfem.VectorPyCoefficient):
        def EvalValue(self, x):
            xi = float(x[0])
            yi = float(x[1])
            zi = 0.0
            if len(x) == 3:
                zi = x[2]
            ret = [- exp(xi)*sin(yi)*cos(zi),
                   - exp(xi)*cos(yi)*cos(zi)]
            if len(x) == 3:
                ret.append(exp(xi)*sin(yi)*sin(zi))
            return ret

    class pFunc_ex(mfem.PyCoefficient):
        def EvalValue(self, x):
            return pFunc_exact(x)

    class fFunc(mfem.VectorPyCoefficient):
        def EvalValue(self, x):
            if len(x) == 3:
                return [0., 0., 0.]
            else:
                return [0., 0.]

    class gFunc(mfem.PyCoefficient):
        def EvalValue(self, x):
            if len(x) == 3:
                return -pFunc_exact(x)
            else:
                return 0.0

    class f_natural(mfem.PyCoefficient):
        def EvalValue(self, x):
            return -pFunc_exact(x)

    device = mfem.Device(device)
    if myid == 0:
        device.Print()

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    if refine == -1:
        ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
        for x in range(ref_levels):
            mesh.UniformRefinement()
    for i in range(refine):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    par_ref_levels = 2
    for l in range(par_ref_levels):
        pmesh.UniformRefinement()

    hdiv_coll = mfem.RT_FECollection(order, dim)
    l2_coll = mfem.L2_FECollection(order, dim)

    R_space = mfem.ParFiniteElementSpace(pmesh, hdiv_coll)
    W_space = mfem.ParFiniteElementSpace(pmesh, l2_coll)

    dimR = R_space.GlobalTrueVSize()
    dimW = W_space.GlobalTrueVSize()

    if verbose:
        print("***********************************************************")
        print("dim(R) = " + str(dimR))
        print("dim(W) = " + str(dimW))
        print("dim(R+W) = " + str(dimR+dimW))
        print("***********************************************************")

    block_offsets = intArray([0, R_space.GetVSize(), W_space.GetVSize()])
    block_offsets.PartialSum()
    block_trueOffsets = intArray([0, R_space.TrueVSize(), W_space.TrueVSize()])
    block_trueOffsets.PartialSum()

    k = mfem.ConstantCoefficient(1.0)

    fcoeff = fFunc(dim)
    fnatcoeff = f_natural()
    gcoeff = gFunc()
    ucoeff = uFunc_ex(dim)
    pcoeff = pFunc_ex()

    x = mfem.BlockVector(block_offsets)

    rhs = mfem.BlockVector(block_offsets)
    trueX = mfem.BlockVector(block_trueOffsets)

    trueRhs = mfem.BlockVector(block_trueOffsets)
    trueRhs.Assign(0.0)

    fform = mfem.ParLinearForm()
    fform.Update(R_space, rhs.GetBlock(0), 0)
    fform.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(fcoeff))
    fform.AddBoundaryIntegrator(
        mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
    fform.Assemble()
    fform.ParallelAssemble(trueRhs.GetBlock(0))

    gform = mfem.ParLinearForm()
    gform.Update(W_space, rhs.GetBlock(1), 0)
    gform.AddDomainIntegrator(mfem.DomainLFIntegrator(gcoeff))
    gform.Assemble()
    gform.ParallelAssemble(trueRhs.GetBlock(1))

    mVarf = mfem.ParBilinearForm(R_space)
    bVarf = mfem.ParMixedBilinearForm(R_space, W_space)

    if pa:
        mVarf.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    mVarf.AddDomainIntegrator(mfem.VectorFEMassIntegrator(k))
    mVarf.Assemble()
    if not pa:
        mVarf.Finalize()

    if pa:
        bVarf.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    bVarf.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator())
    bVarf.Assemble()
    if not pa:
        bVarf.Finalize()

    darcyOp = mfem.BlockOperator(block_trueOffsets)

    if pa:
        opM = mfem.OperatorPtr()
        opB = mfem.OperatorPtr()

        empty_tdof_list = mfem.intArray()
        mVarf.FormSystemMatrix(empty_tdof_list, opM)
        bVarf.FormRectangularSystemMatrix(
            empty_tdof_list, empty_tdof_list, opB)
        Bt = mfem.TransposeOperator(opB.Ptr())

        darcyOp.SetBlock(0, 0, opM.Ptr())
        darcyOp.SetBlock(0, 1, Bt, -1.0)
        darcyOp.SetBlock(1, 0, opB.Ptr(), -1.0)
    else:
        M = mVarf.ParallelAssemble()
        B = bVarf.ParallelAssemble()
        B *= -1
        Bt = mfem.TransposeOperator(B)

        darcyOp.SetBlock(0, 0, M)
        darcyOp.SetBlock(0, 1, Bt)
        darcyOp.SetBlock(1, 0, B)

    if pa:
        Md_PA = mfem.Vector(R_space.GetTrueVSize())
        mVarf.AssembleDiagonal(Md_PA)
        Md_host = mfem.Vector(Md_PA.HostRead(), R_space.GetTrueVSize())
        invMd = mfem.Vector(1.0 / Md_host.GetDataArray())

        BMBt_diag = mfem.Vector(bVarf.Height())
        bVarf.AssembleDiagonal_ADAt(invMd, BMBt_diag)
        ess_tdof_list = mfem.intArray()
        invM = mfem.OperatorJacobiSmoother(Md_PA, ess_tdof_list)
        invS = mfem.OperatorJacobiSmoother(BMBt_diag, ess_tdof_list)
    else:
        Md = mfem.HypreParVector(MPI.COMM_WORLD, M.GetGlobalNumRows(),
                                 M.GetRowStarts())
        M.GetDiag(Md)

        MinvBt = B.Transpose()
        MinvBt.InvScaleRows(Md)
        S = mfem.hypre.ParMult(B, MinvBt)
        invM = mfem.HypreDiagScale(M)
        invS = mfem.HypreBoomerAMG(S)

    invM.iterative_mode = False
    invS.iterative_mode = False

    darcyPr = mfem.BlockDiagonalPreconditioner(block_trueOffsets)
    darcyPr.SetDiagonalBlock(0, invM)
    darcyPr.SetDiagonalBlock(1, invS)

    maxIter = 1000 if pa else 500
    rtol = 1e-6
    atol = 1e-10

    stime = clock()
    solver = mfem.MINRESSolver(MPI.COMM_WORLD)

    solver.SetAbsTol(atol)
    solver.SetRelTol(rtol)
    solver.SetMaxIter(maxIter)
    solver.SetOperator(darcyOp)
    solver.SetPreconditioner(darcyPr)
    solver.SetPrintLevel(1)
    trueX.Assign(0.0)
    solver.Mult(trueRhs, trueX)
    if device.IsEnabled():
        trueX.HostRead()

    solve_time = clock() - stime

    if verbose:
        if solver.GetConverged():
            print("MINRES converged in " + str(solver.GetNumIterations()) +
                  " iterations with a residual norm of " +
                  "{:g}".format(solver.GetFinalNorm()))
        else:
            print("MINRES did not converge in " + str(solver.GetNumIterations()) +
                  " iterations. Residual norm is " +
                  "{:g}".format(solver.GetFinalNorm()))
        print("MINRES solver took " + "{:g}".format(solve_time) + "s.")

    u = mfem.ParGridFunction()
    p = mfem.ParGridFunction()
    u.MakeRef(R_space, x.GetBlock(0), 0)
    p.MakeRef(W_space, x.GetBlock(1), 0)
    u.Distribute(trueX.GetBlock(0))
    p.Distribute(trueX.GetBlock(1))

    order_quad = max(2, 2*order+1)
    irs = [mfem.IntRules.Get(i, order_quad)
           for i in range(mfem.Geometry.NumGeom)]

    err_u = u.ComputeL2Error(ucoeff, irs)
    norm_u = mfem.ComputeGlobalLpNorm(2, ucoeff, pmesh, irs)
    err_p = p.ComputeL2Error(pcoeff, irs)
    norm_p = mfem.ComputeGlobalLpNorm(2, pcoeff, pmesh, irs)

    if verbose:
        print("|| u_h - u_ex || / || u_ex || = " +
              "{:g}".format(err_u / norm_u))
        print("|| p_h - p_ex || / || p_ex || = " +
              "{:g}".format(err_p / norm_p))

    #  Save the refined mesh and the solution in parallel. This output can be
    #  viewed later using GLVis: "glvis -np <np> -m mesh -g sol_*".
    pmesh.Print('mesh.' + smyid, 8)
    u.Save('sol_u.' + smyid, 8)
    p.Save('sol_p.' + smyid, 8)

    # Save data in the Visit format
    visit_dc = mfem.VisItDataCollection("Example5-Parallel", pmesh)
    visit_dc.RegisterField("velocity", u)
    visit_dc.RegisterField("pressure", p)
    visit_format = (mfem.DataCollection.PARALLEL_FORMAT if par_format else
                    mfem.DataCollection.SERIAL_FORMAT)
    visit_dc.SetFormat(visit_format)
    visit_dc.Save()

    # Save data in the ParaVier format
    paraview_dc = mfem.ParaViewDataCollection("Example5P", pmesh)
    paraview_dc.SetPrefixPath("ParaView")
    paraview_dc.SetLevelsOfDetail(order)
    paraview_dc.SetCycle(0)
    paraview_dc.SetDataFormat(mfem.VTKFormat_BINARY)
    paraview_dc.SetHighOrderOutput(True)
    paraview_dc.SetTime(0.0)
    paraview_dc.RegisterField("velocity", u)
    paraview_dc.RegisterField("pressure", p)
    paraview_dc.Save()

    if visualization:
        u_sock = mfem.socketstream("localhost", 19916)
        u_sock << "parallel " << num_procs << " " << myid << "\n"
        u_sock.precision(8)
        u_sock << "solution\n" << pmesh << u << "window_title 'Velocity'\n"
        MPI.Barrier()
        p_sock = mfem.socketstream("localhost", 19916)
        p_sock << "parallel " << num_procs << " " << myid << "\n"
        p_sock.precision(8)
        p_sock << "solution\n" << pmesh << p << "window_title 'Pressure'\n"


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex5 (Darcy Problem)')
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-r', '--refine',
                        action='store', type=int, default=-1,
                        help='Number of times to refine the mesh uniformly')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")
    parser.add_argument('-pf', '--parallel-format',
                        action='store_true',
                        help="Use parallel format for Visit output")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    device = args.device
    pa = args.partial_assembly

    run(order=order,
        meshfile=meshfile,
        visualization=visualization,
        refine=args.refine,
        par_format=args.parallel_format,
        device=device,
        pa=pa)
