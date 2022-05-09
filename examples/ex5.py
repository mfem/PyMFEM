'''
   MFEM example 5

   See c++ version in the MFEM library for more detail 
'''
import os
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp

# time.clock deprecated and removed in PY3.8
import time
try:
    clock = time.process_time
except AttributeError:
    clock = time.clock


def run(order=1,
        meshfile='',
        visualization=False,
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
            return 0.

    class f_natural(mfem.PyCoefficient):
        def EvalValue(self, x):
            return -pFunc_exact(x)

    device = mfem.Device(device)
    device.Print()

    mesh = mfem.Mesh(meshfile, 1, 1)

    dim = mesh.Dimension()

    ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
    for x in range(ref_levels):
        mesh.UniformRefinement()

    hdiv_coll = mfem.RT_FECollection(order, dim)
    l2_coll = mfem.L2_FECollection(order, dim)

    R_space = mfem.FiniteElementSpace(mesh, hdiv_coll)
    W_space = mfem.FiniteElementSpace(mesh, l2_coll)

    dimR = R_space.GetVSize()
    dimW = W_space.GetVSize()

    print("***********************************************************")
    print("dim(R) = " + str(dimR))
    print("dim(W) = " + str(dimW))
    print("dim(R+W) = " + str(dimR+dimW))
    print("***********************************************************")

    block_offsets = intArray([0, dimR, dimW])
    block_offsets.PartialSum()

    k = mfem.ConstantCoefficient(1.0)

    fcoeff = fFunc(dim)
    fnatcoeff = f_natural()
    gcoeff = gFunc()
    ucoeff = uFunc_ex(dim)
    pcoeff = pFunc_ex()

    x = mfem.BlockVector(block_offsets)
    rhs = mfem.BlockVector(block_offsets)

    fform = mfem.LinearForm()
    fform.Update(R_space, rhs.GetBlock(0), 0)
    fform.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(fcoeff))
    fform.AddBoundaryIntegrator(
        mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
    fform.Assemble()

    gform = mfem.LinearForm()
    gform.Update(W_space, rhs.GetBlock(1), 0)
    gform.AddDomainIntegrator(mfem.DomainLFIntegrator(gcoeff))
    gform.Assemble()

    mVarf = mfem.BilinearForm(R_space)
    bVarf = mfem.MixedBilinearForm(R_space, W_space)

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

    darcyOp = mfem.BlockOperator(block_offsets)

    if pa:
        Bt = mfem.TransposeOperator(bVarf)

        darcyOp.SetBlock(0, 0, mVarf)
        darcyOp.SetBlock(0, 1, Bt, -1)
        darcyOp.SetBlock(1, 0, bVarf, -1)
    else:
        M = mVarf.SpMat()
        B = bVarf.SpMat()
        B *= -1
        if mfem.Device.IsEnabled():
            B = mfem.Transpose(B)
        Bt = mfem.TransposeOperator(B)

        darcyOp.SetBlock(0, 0, M)
        darcyOp.SetBlock(0, 1, Bt)
        darcyOp.SetBlock(1, 0, B)

    Md = mfem.Vector(mVarf.Height())
    darcyPrec = mfem.BlockDiagonalPreconditioner(block_offsets)

    if pa:
        mVarf.AssembleDiagonal(Md)
        Md_host = mfem.Vector(Md.HostRead(), mVarf.Height())
        invMd = mfem.Vector(mVarf.Height())

        for i in range(invMd.Size()):
            invMd[i] = 1.0/Md_host[i]

        BMBt_diag = mfem.Vector(bVarf.Height())
        bVarf.AssembleDiagonal_ADAt(invMd, BMBt_diag)
        ess_tdof_list = mfem.intArray()
        invM = mfem.OperatorJacobiSmoother(Md, ess_tdof_list)
        invS = mfem.OperatorJacobiSmoother(BMBt_diag, ess_tdof_list)
    else:
        MinvBt = mfem.Transpose(B)
        M.GetDiag(Md)

        for i in range(Md.Size()):
            MinvBt.ScaleRow(i, 1/Md[i])
        S = mfem.Mult(B, MinvBt)

        invM = mfem.DSmoother(M)
        invS = mfem.GSSmoother(S)

    invM.iterative_mode = False
    invS.iterative_mode = False

    darcyPrec.SetDiagonalBlock(0, invM)
    darcyPrec.SetDiagonalBlock(1, invS)

    maxIter = 1000
    rtol = 1e-6
    atol = 1e-10

    stime = clock()
    solver = mfem.MINRESSolver()
    solver.SetAbsTol(atol)
    solver.SetRelTol(rtol)
    solver.SetMaxIter(maxIter)
    solver.SetOperator(darcyOp)
    solver.SetPreconditioner(darcyPrec)
    solver.SetPrintLevel(1)
    x.Assign(0.0)
    solver.Mult(rhs, x)

    solve_time = clock() - stime

    if solver.GetConverged():
        print("MINRES converged in " + str(solver.GetNumIterations()) +
              " iterations with a residual norm of " + "{:g}".format(solver.GetFinalNorm()))
    else:
        print("MINRES did not converge in " + str(solver.GetNumIterations()) +
              " iterations. Residual norm is " + "{:g}".format(solver.GetFinalNorm()))
    print("MINRES solver took " + str(solve_time) + "s.")

    u = mfem.GridFunction()
    p = mfem.GridFunction()
    u.MakeRef(R_space, x.GetBlock(0), 0)
    p.MakeRef(W_space, x.GetBlock(1), 0)

    order_quad = max(2, 2*order+1)

    irs = [mfem.IntRules.Get(i, order_quad)
           for i in range(mfem.Geometry.NumGeom)]

    norm_p = mfem.ComputeLpNorm(2, pcoeff, mesh, irs)
    norm_u = mfem.ComputeLpNorm(2, ucoeff, mesh, irs)
    err_u = u.ComputeL2Error(ucoeff, irs)
    err_p = p.ComputeL2Error(pcoeff, irs)

    print("|| u_h - u_ex || / || u_ex || = " + "{:g}".format(err_u / norm_u))
    print("|| p_h - p_ex || / || p_ex || = " + "{:g}".format(err_p / norm_p))

    mesh.Print("ex5.mesh", 8)
    u.Save("sol_u.gf", 8)
    p.Save("sol_p.gf", 8)

    visit_dc = mfem.VisItDataCollection("Example5", mesh)
    visit_dc.RegisterField("velocity", u)
    visit_dc.RegisterField("pressure", p)
    visit_dc.Save()

    paraview_dc = mfem.ParaViewDataCollection("Example5", mesh)
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
        sout_u = mfem.socketstream("localhost", 19916)
        sout_u.precision(8)
        sout_u << "solution\n" << mesh << u << "window_title 'Velocity'\n"
        sout_p = mfem.socketstream("localhost", 19916)
        sout_p.precision(8)
        sout_p << "solution\n" << mesh << p << "window_title 'Pressure'\n"


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex5 (Darcy Problem)')
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")

    args = parser.parse_args()
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
        device=device,
        pa=pa)
