'''
   MFEM example 3

   See c++ version in the MFEM library for more detail 
'''
import mfem.ser as mfem
from mfem.ser import intArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, array


def run(order=1,
        static_cond=False,
        freq=1,
        meshfile='',
        visualization=False,
        device='cpu',
        numba=False,
        pa=False):

    kappa = np.pi*freq

    device = mfem.Device(device)
    device.Print()

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    if numba:
        @mfem.jit.vector()
        def E_exact(x, out):
            out[0] = sin(kappa*x[1])
            out[1] = sin(kappa*x[2])
            out[2] = sin(kappa*x[0])

        @mfem.jit.vector()
        def f_exact(x, out):
            out[0] = (1 + kappa**2)*sin(kappa * x[1])
            out[1] = (1 + kappa**2)*sin(kappa * x[2])
            out[2] = (1 + kappa**2)*sin(kappa * x[0])
    else:
        class cE_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                return (sin(kappa * x[1]),
                        sin(kappa * x[2]),
                        sin(kappa * x[0]))
        E_exact = cE_exact()

        class cf_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                return ((1 + kappa**2)*sin(kappa * x[1]),
                        (1 + kappa**2)*sin(kappa * x[2]),
                        (1 + kappa**2)*sin(kappa * x[0]))
        f_exact = cf_exact()

    #   3. Refine the mesh to increase the resolution. In this example we do
    #      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    #      largest number that gives a final mesh with no more than 50,000
    #      elements.

    ref_levels = int(np.floor(np.log(50000./mesh.GetNE())/np.log(2.)/dim))
    for x in range(ref_levels):
        mesh.UniformRefinement()

    #  4. Define a finite element space on the mesh. Here we use the Nedelec
    #     finite elements of the specified order.

    fec = mfem.ND_FECollection(order, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec)

    print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))

    # 5. Determine the list of true (i.e. conforming) essential boundary dofs.
    #    In this example, the boundary conditions are defined by marking all
    #    the boundary attributes from the mesh as essential (Dirichlet) and
    #    converting them to a list of true dofs.

    ess_tdof_list = intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr = intArray(mesh.bdr_attributes.Max())
        ess_bdr = intArray([1]*mesh.bdr_attributes.Max())
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 6. Set up the linear form b(.) which corresponds to the right-hand side
    #    of the FEM linear system, which in this case is (f,phi_i) where f is
    #    given by the function f_exact and phi_i are the basis functions in the
    #    finite element fespace.

    b = mfem.LinearForm(fespace)
    dd = mfem.VectorFEDomainLFIntegrator(f_exact)
    b.AddDomainIntegrator(dd)
    b.Assemble()

    # 7. Define the solution vector x as a finite element grid function
    #    corresponding to fespace. Initialize x by projecting the exact
    #    solution. Note that only values from the boundary edges will be used
    #    when eliminating the non-homogeneous boundary condition to modify the
    #    r.h.s. vector b.

    #from mfem.examples.ex3 import E_exact_cb
    x = mfem.GridFunction(fespace)
    x.ProjectCoefficient(E_exact)

    # 8. Set up the bilinear form corresponding to the EM diffusion operator
    #       curl muinv curl + sigma I, by adding the curl-curl and the mass domain
    #       integrators.

    muinv = mfem.ConstantCoefficient(1.0)
    sigma = mfem.ConstantCoefficient(1.0)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))

    # 9. Assemble the bilinear form and the corresponding linear system,
    #       applying any necessary transformations such as: eliminating boundary
    #       conditions, applying conforming constraints for non-conforming AMR,
    #       static condensation, etc.

    if (static_cond):
        a.EnableStaticCondensation()
    a.Assemble()

    A = mfem.SparseMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
    # Here, original version calls hegith, which is not
    # defined in the header...!?
    print("Size of linear system: " + str(A.Size()))
    import sys
    sys.stdout.flush()
    # 10. Solve
    M = mfem.GSSmoother(A)
    mfem.PCG(A, M, B, X, 1, 500, 1e-12, 0.0)
    sys.stdout.flush()

    # 11. Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, x)

    # 12. Compute and print the L^2 norm of the error.
    import sys
    sys.stdout.write("|| E_h - E ||_{L^2} = " +
                     str(x.ComputeL2Error(E_exact))+"\n")

    mesh.Print('refined.mesh', 8)
    x.Save('sol.gf', 8)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex3 (Definite Maxwell Problem)')
    parser.add_argument('-m', '--mesh',
                        default="beam-tet.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=1.0,
                        help="Set the frequency for the exact")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")
    parser.add_argument('-sc', '--static-condensation',
                        action='store_true',
                        help="Enable static condensation.")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
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
    parser.print_options(args)

    order = args.order
    static_cond = args.static_condensation

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
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
