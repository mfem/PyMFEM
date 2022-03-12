'''
   MFEM example 3 parallel version

   See c++ version in the MFEM library for more detail 
'''
from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from mpi4py import MPI
import numpy as np
from numpy import sin, array


static_cond = False
order = 1

nprc = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)
verbose = (myid == 0)


def run(order=1,
        static_cond=False,
        freq=1,
        meshfile='',
        visualization=False,
        device='cpu',
        numba=False,
        pa=False):

    kappa = np.pi*freq

    # 3. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device)
    if myid == 0:
        device.Print()

    # 4. Read the (serial) mesh from the given mesh file on all processors.  We
    #    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    #    and volume meshes with the same code.
    mesh = mfem.Mesh(meshfile, 1, 1)

    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    if numba:
        @mfem.jit.vector()
        def E_exact(x, out):
            if dim == 3:            
                out[0] = sin(kappa*x[1])
                out[1] = sin(kappa*x[2])
                out[2] = sin(kappa*x[0])
            else:
                out[0] = sin(kappa*x[1])
                out[1] = sin(kappa*x[0])
                if sdim == 3:
                    out[2] = 0.                    
        @mfem.jit.vector()
        def f_exact(x, out):
            if dim == 3:
                out[0] = (1 + kappa**2)*sin(kappa * x[1])
                out[1] = (1 + kappa**2)*sin(kappa * x[2])
                out[2] = (1 + kappa**2)*sin(kappa * x[0])
            else:
                out[0] = (1 + kappa**2)*sin(kappa * x[1])
                out[1] = (1 + kappa**2)*sin(kappa * x[0])
                if sdim == 3:
                    out[2] = 0.
                
    else:
        class cE_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                if dim == 3:
                    return (sin(kappa * x[1]),
                            sin(kappa * x[2]),
                            sin(kappa * x[0]))
                elif sdim == 2:
                        return (sin(kappa * x[1]),
                                sin(kappa * x[0]),)
                else:
                        return (sin(kappa * x[1]),
                                sin(kappa * x[0]),
                                0)
                    
        E_exact = cE_exact()

        class cf_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                if dim == 3:                
                    return ((1 + kappa**2)*sin(kappa * x[1]),
                           (1 + kappa**2)*sin(kappa * x[2]),
                           (1 + kappa**2)*sin(kappa * x[0]))
                elif sdim == 2:
                    return ((1 + kappa**2)*sin(kappa * x[1]),
                           (1 + kappa**2)*sin(kappa * x[0]),)
                else:
                    return ((1 + kappa**2)*sin(kappa * x[1]),
                            (1 + kappa**2)*sin(kappa * x[0]),
                            0.0)
                    
        f_exact = cf_exact()

    # 5. Refine the serial mesh on all processors to increase the resolution. In
    #    this example we do 'ref_levels' of uniform refinement. We choose
    #    'ref_levels' to be the largest number that gives a final mesh with no
    #    more than 1,000 elements

    ref_levels = int(np.floor(np.log(1000./mesh.GetNE())/np.log(2.)/dim))
    for x in range(ref_levels):
        mesh.UniformRefinement()

    # 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution. Once the
    #    parallel mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    par_ref_levels = 2
    for l in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 7. Define a parallel finite element space on the parallel mesh. Here we
    #    use the Nedelec finite elements of the specified order.
    fec = mfem.ND_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    size = fespace.GlobalTrueVSize()

    if verbose:  # note that size should be evaulated on all nodes
        print("Number of finite element unknowns: " + str(size))

    # 8. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, the boundary conditions are defined
    #    by marking all the boundary attributes from the mesh as essential
    #    (Dirichlet) and converting them to a list of true dofs.
    ess_tdof_list = intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr = intArray(mesh.bdr_attributes.Max())
        ess_bdr.Assign(1)
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 9. Set up the parallel linear form b(.) which corresponds to the
    #    right-hand side of the FEM linear system, which in this case is
    #    (f,phi_i) where f is given by the function f_exact and phi_i are the
    #    basis functions in the finite element fespace.
    b = mfem.ParLinearForm(fespace)
    dd = mfem.VectorFEDomainLFIntegrator(f_exact)
    b.AddDomainIntegrator(dd)
    b.Assemble()

    # 10. Define the solution vector x as a parallel finite element grid function
    #     corresponding to fespace. Initialize x by projecting the exact
    #     solution. Note that only values from the boundary edges will be used
    #     when eliminating the non-homogeneous boundary condition to modify the
    #     r.h.s. vector b.
    x = mfem.ParGridFunction(fespace)
    x.ProjectCoefficient(E_exact)

    # 11. Set up the parallel bilinear form corresponding to the EM diffusion
    #     operator curl muinv curl + sigma I, by adding the curl-curl and the
    #     mass domain integrators.
    muinv = mfem.ConstantCoefficient(1.0)
    sigma = mfem.ConstantCoefficient(1.0)
    a = mfem.ParBilinearForm(fespace)
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))

    # 12. Assemble the parallel bilinear form and the corresponding linear
    #     system, applying any necessary transformations such as: parallel
    #     assembly, eliminating boundary conditions, applying conforming
    #     constraints for non-conforming AMR, static condensation, etc
    if (static_cond):
        a.EnableStaticCondensation()
    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

    # 13. Solve the system AX=B using PCG with an AMS preconditioner.
    if pa:
        ams = mfem.MatrixFreeAMS(a, A, fespace, muinv, sigma, None, ess_bdr)
        cg = mfem.CGSolver(MPI.COMM_WORLD)
        cg.SetRelTol(1e-12)
        cg.SetMaxIter(1000)
        cg.SetPrintLevel(1)
        cg.SetOperator(A.Ptr())
        cg.SetPreconditioner(ams)
        cg.Mult(B, X)
    else:
        prec_fespace = (a.SCParFESpace() if a.StaticCondensationIsEnabled()
                        else fespace)
        Amat = A.AsHypreParMatrix()
        
        if verbose:
            print("Size of linear system: " + str(Amat.GetGlobalNumRows()))
        
        ams = mfem.HypreAMS(Amat, prec_fespace)
        pcg = mfem.HyprePCG(Amat)
        pcg.SetTol(1e-12)
        pcg.SetMaxIter(500)
        pcg.SetPrintLevel(2)
        pcg.SetPreconditioner(ams)
        pcg.Mult(B, X)

    # 14. Recover the parallel grid function corresponding to X. This is the
    #     local finite element solution on each processor.
    a.RecoverFEMSolution(X, b, x)

    # 15. Compute and print the L^2 norm of the error.
    err = x.ComputeL2Error(E_exact)
    if verbose:  # note that err should be evaulated on all nodes
        print("|| E_h - E ||_{L^2} = " + "{:g}".format(err))

    # 16. Save the refined mesh and the solution in parallel. This output can
    #     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    x.Save('sol.'+smyid)
    pmesh.Print('mesh.'+smyid)

    # 17. Send the solution by socket to a GLVis server
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "parallel " << nprc << " " << myid << "\n"
        sol_sock << "solution\n" << pmesh << x
        sol_sock << "window_title 'Soluiton real part'"
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex3 (Definite Maxwell Problem)')
    parser.add_argument('-m', '--mesh',
                        default="beam-tet.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=1.0,
                        help="Set the frequency for the exact")
    parser.add_argument('-sc', '--static-condensation',
                        action='store_true',
                        help="Enable static condensation.")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    try:
        from numba import jit
        HAS_NUMBA = True
    except ImportError:
        HAS_NUMBA = False
    parser.add_argument("-n", "--numba",
                        default=int(HAS_NUMBA),
                        type=int,
                        help="Use Number compiled coefficient")

    args = parser.parse_args()
    args.numba = bool(args.numba)

    if verbose:
        parser.print_options(args)

    order = args.order
    static_cond = args.static_condensation

    from os.path import expanduser, join, dirname
    meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
    
    visualization = args.visualization
    device = args.device
    pa = args.partial_assembly
    freq = args.frequency
    numba = args.numba

    run(freq=freq,
        order=order,
        static_cond=static_cond,
        meshfile=meshfile,
        visualization=visualization,
        device=device,
        pa=pa,
        numba=numba)
