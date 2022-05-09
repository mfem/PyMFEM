
'''
   MFEM example 24
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array, floor, log

freq = 1.0


def run(order=1,
        prob=0,
        static_cond=False,
        mesh_file="",
        visualization=1,
        device="cpu",
        pa=False,
        numba=True):

    kappa = freq * pi

    # 2. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device)
    device.Print()

    # 3. Read the mesh from the given mesh file. We can handle triangular,
    #    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
    #    the same code.
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    #    largest number that gives a final mesh with no more than 50,000
    #    elements.
    ref_levels = int(floor(log(50000./mesh.GetNE())/log(2.)/dim))
    for l in range(ref_levels):
        mesh.UniformRefinement()
    mesh.ReorientTetMesh()

    # 5. Define a finite element space on the mesh. Here we use Nedelec or
    #   Raviart-Thomas finite elements of the specified order.
    if prob == 0:
        trial_fec = mfem.H1_FECollection(order, dim)
        test_fec = mfem.ND_FECollection(order, dim)
    elif prob == 1:
        trial_fec = mfem.ND_FECollection(order, dim)
        test_fec = mfem.RT_FECollection(order-1, dim)
    else:
        trial_fec = mfem.RT_FECollection(order-1, dim)
        test_fec = mfem.L2_FECollection(order-1, dim)

    trial_fes = mfem.FiniteElementSpace(mesh, trial_fec)
    test_fes = mfem.FiniteElementSpace(mesh, test_fec)

    trial_size = trial_fes.GetTrueVSize()
    test_size = test_fes.GetTrueVSize()

    if prob == 0:
        print("Number of Nedelec finite element unknowns: " + str(test_size))
        print("Number of H1 finite element unknowns: " + str(trial_size))
    elif prob == 1:
        print("Number of Nedelec finite element unknowns: " + str(test_size))
        print("Number of Raviart-Thomas finite element unknowns: " + str(trial_size))
    else:
        print("Number of Raviart-Thomas finite element unknowns: " + str(test_size))
        print("Number of L2 finite element unknowns: " + str(trial_size))

    # 6. Define the solution vector as a finite element grid function
    #    corresponding to the trial fespace.
    gftest = mfem.GridFunction(test_fes)
    gftrial = mfem.GridFunction(trial_fes)
    x = mfem.GridFunction(test_fes)

    if numba:
        @mfem.jit.scalar()
        def p_coef(x):
            if dim == 3:
                return sin(x[0]) * sin(x[1]) * sin(x[2])
            elif dim == 2:
                return sin(x[0]) * sin(x[1])
            return 0.0

        @mfem.jit.vector(sdim=sdim, vdim=dim)
        def gradp_coef(x, out, sdim_, dim_):
            if sdim_ == 3:
                out[0] = cos(x[0]) * sin(x[1]) * sin(x[2])
                out[1] = sin(x[0]) * cos(x[1]) * sin(x[2])
                out[2] = sin(x[0]) * sin(x[1]) * cos(x[2])
            elif sdim_ == 2:
                out[0] = cos(x[0]) * sin(x[1])
                out[1] = sin(x[0]) * cos(x[1])
                if dim_ == 3:
                    out[2] = 0.0

        @mfem.jit.vector(sdim=sdim, vdim=dim)
        def v_coef(x, out, sdim_, dim_):
            if sdim_ == 3:
                out[0] = sin(kappa * x[1])
                out[1] = sin(kappa * x[2])
                out[2] = sin(kappa * x[0])
            else:
                out[0] = sin(kappa * x[1])
                out[1] = sin(kappa * x[0])
                if dim_ == 3:
                    out[2] = 0.0

        @mfem.jit.scalar()
        def cdiv_gradp_coeff(x):
            if dim == 3:
                return -3.0 * sin(x[0]) * sin(x[1]) * sin(x[2])
            elif dim == 2:
                return -2.0 * sin(x[0]) * sin(x[1])
            return 0.0

        @mfem.jit.vector(sdim=sdim, vdim=dim)
        def curlv_coef(x, out, sdim_, dim_):
            if sdim_ == 3:
                out[0] = -kappa * cos(kappa * x[2])
                out[1] = -kappa * cos(kappa * x[0])
                out[2] = -kappa * cos(kappa * x[1])
            else:
                for i in range(dim_):
                    out[i] = 0.0

        @mfem.jit.scalar()
        def divgradp_coef(x):
            if dim == 3:
                return -3.0 * sin(x[0]) * sin(x[1]) * sin(x[2])
            elif dim == 2:
                return -2.0 * sin(x[0]) * sin(x[1])
            return 0.0

    else:
        class cp_exact(mfem.PyCoefficient):
            def EvalValue(self, x):
                if dim == 3:
                    return sin(x[0]) * sin(x[1]) * sin(x[2])
                elif dim == 2:
                    return sin(x[0]) * sin(x[1])
                return 0.0

        class cgradp_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                if dim == 3:
                    return [cos(x[0]) * sin(x[1]) * sin(x[2]),
                            sin(x[0]) * cos(x[1]) * sin(x[2]),
                            sin(x[0]) * sin(x[1]) * cos(x[2]), ]
                elif dim == 2:
                    v1 = cos(x[0]) * sin(x[1])
                    v2 = sin(x[0]) * cos(x[1])
                    if len(x) == 3:
                        return v1, v2, 0.0
                    else:
                        return v1, v2

        class cdiv_gradp_exact(mfem.PyCoefficient):
            def EvalValue(self, x):
                if dim == 3:
                    return -3.0 * sin(x[0]) * sin(x[1]) * sin(x[2])
                elif dim == 2:
                    return -2.0 * sin(x[0]) * sin(x[1])
                return 0.0

        class cv_exact(mfem.VectorPyCoefficient):
            def __init__(self):
                mfem.VectorPyCoefficient.__init__(self, sdim)

            def EvalValue(self, x):
                if dim == 3:
                    return [sin(kappa * x[1]),
                            sin(kappa * x[2]),
                            sin(kappa * x[0])]
                else:
                    v1 = sin(kappa * x[1])
                    v2 = sin(kappa * x[0])
                    if len(x) == 3:
                        return (v1, v2, 0.0)
                    else:
                        return (v1, v2)

        class ccurlv_exact(mfem.PyCoefficient):
            def EvalValue(self, x):
                if dim == 3:
                    return [-kappa * cos(kappa * x[2]),
                            -kappa * cos(kappa * x[0]),
                            -kappa * cos(kappa * x[1]), ]
                else:
                    return [0.0] * self.vdim

        p_coef = cp_exact()
        gradp_coef = cgradp_exact()
        v_coef = cv_exact()
        curlv_coef = ccurlv_exact()
        divgradp_coef = cdiv_gradp_exact()

    if prob == 0:
        gftrial.ProjectCoefficient(p_coef)
    elif prob == 1:
        gftrial.ProjectCoefficient(v_coef)
    else:
        gftrial.ProjectCoefficient(gradp_coef)

    gftrial.SetTrueVector()
    gftrial.SetFromTrueVector()

    # 7. Set up the bilinear forms for L2 projection.
    one = mfem.ConstantCoefficient(1.0)
    a = mfem.BilinearForm(test_fes)
    a_mixed = mfem.MixedBilinearForm(trial_fes, test_fes)
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
        a_mixed.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

    if prob == 0:
        a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(one))
        a_mixed.AddDomainIntegrator(mfem.MixedVectorGradientIntegrator(one))
    elif prob == 1:
        a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(one))
        a_mixed.AddDomainIntegrator(mfem.MixedVectorCurlIntegrator(one))
    else:
        a.AddDomainIntegrator(mfem.MassIntegrator(one))
        a_mixed.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator(one))

    # 8. Assemble the bilinear form and the corresponding linear system,
    #    applying any necessary transformations such as: eliminating boundary
    #    conditions, applying conforming constraints for non-conforming AMR,
    #    static condensation, etc.
    if static_cond:
        a.EnableStaticCondensation()

    a.Assemble()
    if not pa:
        a.Finalize()

    a_mixed.Assemble()
    if not pa:
        a_mixed.Finalize()

    if pa:
        a_mixed.Mult(gftrial, x)
    else:
        mixed = a_mixed.SpMat()
        mixed.Mult(gftrial, x)

    # 9. Define and apply a PCG solver for Ax = b with Jacobi preconditioner.
    rhs = mfem.GridFunction(test_fes)
    rhs.Assign(x)
    x.Assign(0.0)

    cg = mfem.CGSolver()
    cg.SetRelTol(1e-12)
    cg.SetMaxIter(1000)
    cg.SetPrintLevel(1)
    if pa:
        ess_tdof_list = mfem.intArray()
        Jacobi = mfem.OperatorJacobiSmoother(a, ess_tdof_list)
        cg.SetOperator(a)
        cg.SetPreconditioner(Jacobi)
        cg.Mult(rhs, x)
    else:
        Amat = a.SpMat()
        Jacobi = mfem.DSmoother(Amat)
        cg.SetOperator(Amat)
        cg.SetPreconditioner(Jacobi)
        cg.Mult(rhs, x)

    # 10. Compute the same field by applying a DiscreteInterpolator.
    discreteInterpolant = mfem.GridFunction(test_fes)
    dlo = mfem.DiscreteLinearOperator(trial_fes, test_fes)
    if prob == 0:
        dlo.AddDomainInterpolator(mfem.GradientInterpolator())
    elif prob == 1:
        dlo.AddDomainInterpolator(mfem.CurlInterpolator())
    else:
        dlo.AddDomainInterpolator(mfem.DivergenceInterpolator())

    dlo.Assemble()
    dlo.Mult(gftrial, discreteInterpolant)

    # 11. Compute the projection of the exact field.
    exact_proj = mfem.GridFunction(test_fes)
    if prob == 0:
        exact_proj.ProjectCoefficient(gradp_coef)
    elif prob == 1:
        exact_proj.ProjectCoefficient(curlv_coef)
    else:
        exact_proj.ProjectCoefficient(divgradp_coef)

    exact_proj.SetTrueVector()
    exact_proj.SetFromTrueVector()

    print()  # put empty line
    # 12. Compute and print the L_2 norm of the error.
    if prob == 0:
        errSol = x.ComputeL2Error(gradp_coef)
        errInterp = discreteInterpolant.ComputeL2Error(gradp_coef)
        errProj = exact_proj.ComputeL2Error(gradp_coef)

        print(" Solution of (E_h,v) = (grad p_h,v) for E_h and v in H(curl): " +
              "|| E_h - grad p ||_{L_2} = " + "{:g}".format(errSol) + "\n")
        print(" Gradient interpolant E_h = grad p_h in H(curl): || E_h - grad p" +
              " ||_{L_2} = " + "{:g}".format(errInterp) + "\n")
        print(" Projection E_h of exact grad p in H(curl): || E_h - grad p " +
              "||_{L_2} = " + "{:g}".format(errProj) + "\n")
    elif prob == 1:
        errSol = x.ComputeL2Error(curlv_coef)
        errInterp = discreteInterpolant.ComputeL2Error(curlv_coef)
        errProj = exact_proj.ComputeL2Error(curlv_coef)

        print(" Solution of (E_h,w) = (curl v_h,w) for E_h and w in H(div): " +
              "|| E_h - curl v ||_{L_2} = " + "{:g}".format(errSol) + "\n")
        print(" Curl interpolant E_h = curl v_h in H(div): || E_h - curl v " +
              "||_{L_2} = " + "{:g}".format(errInterp) + "\n")
        print(" Projection E_h of exact curl v in H(div): || E_h - curl v " +
              "||_{L_2} = " + "{:g}".format(errProj) + "\n")
    else:
        order_quad = max(2, 2*order+1)
        irs = [mfem.IntRules.Get(i, order_quad)
               for i in range(mfem.Geometry.NumGeom)]

        errSol = x.ComputeL2Error(divgradp_coef, irs)
        errInterp = discreteInterpolant.ComputeL2Error(divgradp_coef, irs)
        errProj = exact_proj.ComputeL2Error(divgradp_coef, irs)

        print(" Solution of (f_h,q) = (div v_h,q) for f_h and q in L_2: " +
              "|| f_h - div v ||_{L_2} = " + "{:g}".format(errSol) + "\n")
        print(" Divergence interpolant f_h = div v_h in L_2: || f_h - div v " +
              "||_{L_2} = " + "{:g}".format(errInterp) + "\n")
        print(" Projection f_h of exact div v in L_2: || f_h - div v " +
              "||_{L_2} = " + "{:g}".format(errProj) + "\n")

    # 13. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g sol.gf".
    mesh.Print("refined.mesh", 8)
    x.Save("sol.gf", 8)

    # 14. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << x
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex24 (Mixed finite element spaces)')

    parser.add_argument('-m', '--mesh',
                        default='beam-hex.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-p', '--problem-type',
                        action='store', type=int, default=0,
                        help="Choose between 0: grad, 1: curl, 2: div")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
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

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    numba = (args.numba == 1)

    run(order=args.order,
        prob=args.problem_type,
        static_cond=args.static_condensation,
        mesh_file=meshfile,
        visualization=args.visualization,
        device=args.device,
        pa=args.partial_assembly,
        numba=numba)
