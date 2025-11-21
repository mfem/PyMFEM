'''
   MFEM example 5

   See c++ version in the MFEM library for more detail
'''
import os
import time
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp

def run(order=1,
        meshfile='',
        nx=0,
        ny=0,
        dg=False,
        brt=False,
        td=False,
        hb=False,
        rd=False,
        sn=False,
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

    if meshfile != "":
        mesh = mfem.Mesh(meshfile, 1, 1)
    else:
        if ny <= 0:
            ny = nx
        mesh = mfem.Mesh.MakeCartesian2D(nx, ny, mfem.Element.QUADRILATERL)

    dim = mesh.Dimension()

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    #    largest number that gives a final mesh with no more than 10,000
    #    elements.
    
    ref_levels = int(np.floor(np.log(10000./mesh.GetNE())/np.log(2.)/dim))
    for x in range(ref_levels):
        mesh.UniformRefinement()

    # 5. Define a finite element space on the mesh. Here we use the
    #    Raviart-Thomas finite elements of the specified order.
    if dg:
        R_coll = mfem.L2_FECollection(order, dim, mfem.BasisType.GaussLobatto)
    elif brt:
        R_coll = mfem.BrokenRT_FECollection(order, dim)
        R_coll_dg = mfem.L2_FECollection(order+1, dim)
    else:
        R_coll = mfem.RT_FECollection(order, dim)
    
    W_coll = mfem.L2_FECollection(order, dim)
    
    R_space = mfem.FiniteElementSpace(mesh, R_coll, dim if dg else 1)
    if brt:
        R_space_hg = mfem.FiniteElementSpace(mesh, R_coll_dg, dim)
    W_space = mfem.FiniteElementSpace(mesh, W_coll)

    darcy = mfem.DarcyForm(R_space, W_space)

    # 6. Define the BlockStructure of the problem, i.e. define the array of
    #    offsets for each variable. The last component of the Array is the sum
    #    of the dimensions of each block.

    block_offsets = darcy.GetOffsets().ToList()

    print("***********************************************************")
    print("dim(R) = "  + str(block_offsets[1] - block_offsets[0]))
    print("dim(W) = "  + str(block_offsets[2] - block_offsets[1]))
    print("dim(R+W) = "  + str(block_offsets[-1]))
    print("***********************************************************")

    # 7. Define the coefficients, analytical solution, and rhs of the PDE.

    kcoeff = mfem.ConstantCoefficient(1.0) # acoustic resistance
    ikcoeff = mfem.RatioCoefficient(1., kcoeff) # inverse acoustic resistance

    fcoeff = fFunc(dim)
    fnatcoeff = f_natural()
    gcoeff = gFunc()
    ucoeff = uFunc_ex(dim)
    pcoeff = pFunc_ex()

    # 8. Allocate memory (x, rhs) for the analytical solution and the right hand
    #    side.  Define the GridFunction u,p for the finite element solution and
    #    linear forms fform and gform for the right hand side.  The data
    #    allocated by x and rhs are passed as a reference to the grid functions
    #    (u,p) and the linear forms (fform, gform).
    
    mt = device.GetMemoryType();
    x = mfem.BlockVector(block_offsets, mt)
    rhs = mfem.BlockVector(block_offsets, mt)
    
    fform = mfem.LinearForm()
    fform.Update(R_space, rhs.GetBlock(0), 0)

    if dg:
        fform.AddDomainIntegrator(mfem.VectorDomainLFIntegrator(fcoeff))
        fform.AddBdrFaceIntegrator(mfem.VectorBoundaryFluxLFIntegrator(fnatcoeff))
   else
   {
      fform.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(fcoeff))
      if brt:
         fform.AddBdrFaceIntegrator(mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
      else:
         fform.AddBoundaryIntegrator(mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
    

    fform.Assemble()
    fform.SyncAliasMemory(rhs)
       
    gform = mfem.LinearForm()
    gform.Update(W_space, rhs.GetBlock(1), 0)
    gform.AddDomainIntegrator(mfem.DomainLFIntegrator(gcoeff))
    gform.Assemble()
    gform.SyncAliasMemory(rhs)
       
    # 9. Assemble the finite element matrices for the Darcy operator
    #                            D = [ M  B^T ]
    #                                [ B   0  ]
    #     where:
    #
    #     M = \int_\Omega k u_h \cdot v_h d\Omega   u_h, v_h \in R_h
    #     B   = -\int_\Omega \div u_h q_h d\Omega   u_h \in R_h, q_h \in W_h
    mVarf = darcy.GetFluxMassForm()
    bVarf = darcy.GetFluxDivForm()

    if dg:
         mtVarf = darcy.GetPotentialMassForm()        
         mVarf.AddDomainIntegrator(mfem.VectorMassIntegrator(kcoeff))
         bVarf.AddDomainIntegrator(mfem.VectorDivergenceIntegrator())
         bVarf.AddInteriorFaceIntegrator(mfem.TransposeIntegrator(
                                          mfem.DGNormalTraceIntegrator(-1.)))
         mtVarf.AddInteriorFaceIntegrator(mfem.HDGDiffusionIntegrator(ikcoeff, td))
    else:
        mVarf.AddDomainIntegrator(mfem.VectorFEMassIntegrator(kcoeff))
        bVarf.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator)
       if brt:
           bVarf.AddInteriorFaceIntegrator(mfem.TransposeIntegrator(
                                             mfem.DGNormalTraceIntegrator(-1.)))

    #set hybridization / assembly level
    chrono = time.time()
       
    ess_flux_tdofs_list = mfem.intArray()
    if hybridization:
       trace_coll = mfem.DG_Interface_FECollection(order, dim);
       trace_space = mfem.FiniteElementSpace(mesh, trace_coll);
       darcy.EnableHybridization(trace_space,
                                 mfem.NormalTraceJumpIntegrator(),
                                 ess_flux_tdofs_list)
    elif (reduction and (dg or brt)):
       darcy.EnableFluxReduction()
       
    if pa:
        mVarf.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

    darcy.Assemble()

    pDarcyOp = mfem.OperatorHandle()
    X = mfem.Vector()
    B = mfem.Vector()       
    x.Assign(0.0)
    darcy.FormLinearSystem(ess_flux_tdofs_list, x, rhs, pDarcyOp, X, B)

    print("Assembly took " + str(time.time() - chrono))
       
    maxIter = 1000
    rtol = 1e-6
    atol = 1e-10

   if (hybridization or (reduction and (dg or brt))):
       # 10. Construct the preconditioner       
       prec = mfem.GSSmoother(pDarcyOp.AsSparseMatrix())
       
       # 11. Solve the linear system with GMRES.
       #     Check the norm of the unpreconditioned residual.
       chrono = time.time()
       solver = mfem.GMRESSolver()
       solver.SetAbsTol(atol)
       solver.SetRelTol(rtol)
       solver.SetMaxIter(maxIter)
       solver.SetOperator(pDarcyOp)
       solver.SetPreconditioner(prec)
       solver.SetPrintLevel(1)

       p


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
                        default="",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-nx", "--ncells-x",
                        action='store', type=int, default=0,
                        help="Number of cells in x.")
    parser.add_argument("-ny", "--ncells-y",
                        action='store', type=int, default=0,
                        help="Number of cells in y.")
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).");
    parser.add_argument("-dg", "--discontinuous",
                        action='store_true',
                        help="Enable DG elements for fluxes.");
    parser.add_argument("-brt", "--broken-RT",
                        action='store_true',
                        help="Enable broken RT elements for fluxes.")
    parser.add_argument("-td", "--stab_diff",
                        action='store', type=float, default=0.5,
                        help="Diffusion stabilization factor (1/2=default)")
    parser.add_argument("-hb", "--hybridization",
                        action='store_true',
                        help="Enable hybridization.")
    parser.add_argument("-rd", "--reduction",
                        action='store_true',
                        help="Enable reduction of DG flux.")
    parser.add_argument("-sn", "--solution-norm",
                        action='store_true',
                        help="Solution norm (0=native, 1=flux, 2=potential).")
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
    parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    device = args.device
    pa = args.partial_assembly

    run(order=order,
        meshfile=meshfile,
        nx=args.ncells_x,
        ny=args.ncells_y,        
        dg=args.discontinous,
        brt=args.broken_RF,
        td=args.stab_diff,
        hb=args.hybridizatoin,
        rd=args.reduciton,
        sn=args.solution_norm,
        visualization=visualization,
        device=device,
        pa=pa)
