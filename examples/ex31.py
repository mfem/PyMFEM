'''
   MFEM example 31

   See c++ version in the MFEM library for more detail

   Sample runs:  python ex31.py -m ../data/inline-segment.mesh -o 2
                 python ex31.py -m ../data/hexagon.mesh -o 2
                 python ex31.py -m ../data/star.mesh -o 2
                 python ex31.py -m ../data/fichera.mesh -o 3 -r 1
                 python ex31.py -m ../data/square-disc-nurbs.mesh -o 3
                 python ex31.py -m ../data/amr-quad.mesh -o 2 -r 1
                 python ex31.py -m ../data/amr-hex.mesh -r 1
'''
import mfem.ser as mfem
from mfem.ser import intArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt

sqrt1_2 = 1/sqrt(2)
sqrt2 = sqrt(2)


def run(order=1,
        refine=2,
        freq=1,
        meshfile='',
        visualization=False,
        numba=False):

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    # 3. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement (2 by default, or specified on
    #    the command line with -r)
    for x in range(refine):
        mesh.UniformRefinement()

    # 4. Define a finite element space on the mesh. Here we use the Nedelec
    #    finite elements of the specified order restricted to 1D, 2D, or 3D
    #    depending on the dimension of the given mesh file.
    if dim == 1:
        fec = mfem.ND_R1D_FECollection(order, dim)
    elif dim == 2:
        fec = mfem.ND_R2D_FECollection(order, dim)
    else:
        fec = mfem.ND_FECollection(order, dim)

    fespace = mfem.FiniteElementSpace(mesh, fec)
    print("Number of H(curl) unknowns: " + str(fespace.GetTrueVSize()))

    # 5. Determine the list of true essential boundary dofs. In this example,
    #    the boundary conditions are defined by marking all the boundary
    #    attributes from the mesh as essential (Dirichlet) and converting them
    #    to a list of true dofs.
    ess_tdof_list = intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr = intArray([1]*mesh.bdr_attributes.Max())
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 5.b Define coefficent to use. Here, we define outside this function.
    #     We call decoratr function (mfem.jit.xxx) manually. Also, note
    #     we are passing dim as params, so that numba uses proper dim
    #     when compiling it.
    if numba:
        params = {"dim": dim, "kappa": pi*freq}
        E = mfem.jit.vector(sdim=3, params=params)(E_exact)
        CurlE = mfem.jit.vector(sdim=3, params=params)(CurlE_exact)
        f = mfem.jit.vector(sdim=3, params=params)(f_exact)
    else:
        pass  # ToDo provide regular example.

    # 6. Set up the linear form b(.) which corresponds to the right-hand side
    #    of the FEM linear system, which in this case is (f,phi_i) where f is
    #    given by the function f_exact and phi_i are the basis functions in
    #    the finite element fespace.

    b = mfem.LinearForm(fespace)
    b.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(f))
    b.Assemble()

    # 7. Define the solution vector x as a finite element grid function
    #    corresponding to fespace. Initialize x by projecting the exact
    #    solution. Note that only values from the boundary edges will be used
    #    when eliminating the non-homogeneous boundary condition to modify the
    #    r.h.s. vector b.
    sol = mfem.GridFunction(fespace)
    sol.ProjectCoefficient(E)

    # 8. Set up the bilinear form corresponding to the EM diffusion
    #    operator curl muinv curl + sigma I, by adding the curl-curl and the
    #    mass domain integrators.
    mat = array([[2.0, sqrt1_2, 0, ],
                 [sqrt1_2, 2.0, sqrt1_2],
                 [0.0, sqrt1_2, 2.0, ], ])
    sigmaMat = mfem.DenseMatrix(mat)

    muinv = mfem.ConstantCoefficient(1.0)
    sigma = mfem. MatrixConstantCoefficient(sigmaMat)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))

    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, sol, b, A, X, B)

    # 10. Solve the system A X = B.
    AM = A.AsSparseMatrix()
    M = mfem.GSSmoother(AM)
    mfem.PCG(A, M, B, X, 1, 500, 1e-12, 0.0)

    # 12. Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, sol)

    # 13. Compute and print the CurlE norm of the error.
    print("|| E_h - E ||_{Hcurl} = " +
          "{:g}".format(sol.ComputeHCurlError(E, CurlE))+"\n")

    # 14. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g sol.gf"
    mesh.Print('refined.mesh', 8)
    sol.Save('sol.gf', 8)

    # 15. Send the solution by socket to a GLVis server.
    if visualization:

        solCoef = mfem.VectorGridFunctionCoefficient(sol)
        dsolCoef = mfem.CurlGridFunctionCoefficient(sol)

        if dim == 1:
            x_sock = mfem.socketstream("localhost", 19916)
            y_sock = mfem.socketstream("localhost", 19916)
            z_sock = mfem.socketstream("localhost", 19916)
            dy_sock = mfem.socketstream("localhost", 19916)
            dz_sock = mfem.socketstream("localhost", 19916)

            x_sock.precision(8)
            y_sock.precision(8)
            z_sock.precision(8)
            dy_sock.precision(8)
            dz_sock.precision(8)

            xVec = mfem.Vector([1, 0, 0.])
            yVec = mfem.Vector([0, 1, 0.])
            zVec = mfem.Vector([0, 0, 1.])

            xVecCoef = mfem.VectorConstantCoefficient(xVec)
            yVecCoef = mfem.VectorConstantCoefficient(yVec)
            zVecCoef = mfem.VectorConstantCoefficient(zVec)

            fec_h1 = mfem.H1_FECollection(order, dim)
            fec_l2 = mfem.L2_FECollection(order-1, dim)

            fes_h1 = mfem.FiniteElementSpace(mesh, fec_h1)
            fes_l2 = mfem.FiniteElementSpace(mesh, fec_l2)

            xComp = mfem.GridFunction(fes_l2)
            yComp = mfem.GridFunction(fes_h1)
            zComp = mfem.GridFunction(fes_h1)

            dyComp = mfem.GridFunction(fes_l2)
            dzComp = mfem.GridFunction(fes_l2)

            xCoef = mfem.InnerProductCoefficient(xVecCoef, solCoef)
            yCoef = mfem.InnerProductCoefficient(yVecCoef, solCoef)
            zCoef = mfem.InnerProductCoefficient(zVecCoef, solCoef)

            xComp.ProjectCoefficient(xCoef)
            yComp.ProjectCoefficient(yCoef)
            zComp.ProjectCoefficient(zCoef)

            x_sock << "solution\n" << mesh << xComp
            x_sock.flush()
            x_sock << "window_title 'X component'"
            x_sock.flush()
            y_sock << "solution\n" << mesh << yComp
            y_sock.flush()
            y_sock << "window_geometry 403 0 400 350 " << "window_title 'Y component'"
            y_sock.flush()
            z_sock << "solution\n" << mesh << zComp
            z_sock.flush()
            z_sock << "window_geometry 806 0 400 350 " << "window_title 'Z component'"
            z_sock.flush()

            dyCoef = mfem.InnerProductCoefficient(yVecCoef, dsolCoef)
            dzCoef = mfem.InnerProductCoefficient(zVecCoef, dsolCoef)
            dyComp.ProjectCoefficient(dyCoef)
            dzComp.ProjectCoefficient(dzCoef)

            dy_sock << "solution\n" << mesh << dyComp
            dy_sock.flush()
            dy_sock << "window_geometry 403 375 400 350 " << "window_title 'Y component of Curl'"
            dy_sock.flush()

            dy_sock << "solution\n" << mesh << dzComp
            dy_sock.flush()
            dy_sock << "window_geometry 403 375 400 350 " << "window_title 'Z component of Curl'"
            dy_sock.flush()

        elif dim == 2:
            xy_sock = mfem.socketstream("localhost", 19916)
            z_sock = mfem.socketstream("localhost", 19916)
            dxy_sock = mfem.socketstream("localhost", 19916)
            dz_sock = mfem.socketstream("localhost", 19916)

            mat = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
            xyMat = mfem.DenseMatrix(mat)
            xyMatCoef = mfem.MatrixConstantCoefficient(xyMat)
            zVec = array([0, 0, 1])
            zVecCoef = mfem.VectorConstantCoefficient(zVec)

            xyCoef = mfem.MatrixVectorProductCoefficient(xyMatCoef, solCoef)
            zCoef = mfem.InnerProductCoefficient(zVecCoef, solCoef)

            fec_h1 = mfem.H1_FECollection(order, dim)
            fec_nd = mfem.ND_FECollection(order, dim)
            fec_rt = mfem.RT_FECollection(order-1, dim)
            fec_l2 = mfem.L2_FECollection(order-1, dim)

            fes_h1 = mfem.FiniteElementSpace(mesh, fec_h1)
            fes_nd = mfem.FiniteElementSpace(mesh, fec_nd)
            fes_rt = mfem.FiniteElementSpace(mesh, fec_rt)
            fes_l2 = mfem.FiniteElementSpace(mesh, fec_l2)

            xyComp = mfem.GridFunction(fes_nd)
            zComp = mfem.GridFunction(fes_h1)

            dxyComp = mfem.GridFunction(fes_rt)
            dzComp = mfem.GridFunction(fes_l2)

            xyComp.ProjectCoefficient(xyCoef)
            zComp.ProjectCoefficient(zCoef)

            xy_sock.precision(8)
            xy_sock << "solution\n" << mesh << xyComp
            xy_sock << "window_title 'XY components'\n"
            xy_sock.flush()
            z_sock << "solution\n" << mesh << zComp
            z_sock << "window_geometry 403 0 400 350 " << "window_title 'Z component'"
            z_sock.flush()

            dxyCoef = mfem.MatrixVectorProductCoefficient(xyMatCoef, dsolCoef)
            dzCoef = mfem.InnerProductCoefficient(zVecCoef, dsolCoef)

            dxyComp.ProjectCoefficient(dxyCoef)
            dzComp.ProjectCoefficient(dzCoef)

            dxy_sock << "solution\n" << mesh << dxyComp
            dxy_sock << "window_geometry 0 375 400 350 " << "window_title 'XY components of Curl'"
            dxy_sock.flush()
            dz_sock << "solution\n" << mesh << dzComp
            dz_sock << "window_geometry 403 375 400 350 " << "window_title 'Z component of Curl'"
            dz_sock.flush()

        else:
            sol_sock = mfem.socketstream("localhost", 19916)
            dsol_sock = mfem.socketstream("localhost", 19916)

            fec_rt = mfem.RT_FECollection(order-1, dim)
            fes_rt = mfem.FiniteElementSpace(mesh, fec_rt)

            dsol = mfem.GridFunction(fes_rt)

            dsol.ProjectCoefficient(dsolCoef)

            sol_sock.precision(8)
            dsol_sock.precision(8)
            sol_sock << "solution\n" << mesh << sol << "window_title 'Solution'"
            sol_sock.flush()
            dsol_sock << "solution\n" << mesh << dsol
            dsol_sock.flush()
            dsol_sock << "window_geometry 0 375 400 350 " << "window_title 'Curl of solution'"
            dsol_sock.flush()


def E_exact(x, E):
    if dim == 1:
        E[0] = 1.1 * sin(kappa * x[0] + 0.0 * pi)
        E[1] = 1.2 * sin(kappa * x[0] + 0.4 * pi)
        E[2] = 1.3 * sin(kappa * x[0] + 0.9 * pi)
    elif dim == 2:
        E[0] = 1.1 * sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        E[1] = 1.2 * sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        E[2] = 1.3 * sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)
    else:
        E[0] = 1.1 * sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        E[1] = 1.2 * sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        E[2] = 1.3 * sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)

        for i in range(3):
            E[i] = E[i]*cos(kappa * x[2])


def CurlE_exact(x, dE):
    if dim == 1:
        c4 = cos(kappa * x[0] + 0.4 * pi)
        c9 = cos(kappa * x[0] + 0.9 * pi)

        dE[0] = 0.0
        dE[1] = -1.3 * c9
        dE[2] = 1.2 * c4
        for i in range(3):
            dE[i] = dE[i]*kappa

    elif dim == 2:
        c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)

        dE[0] = 1.3 * c9
        dE[1] = -1.3 * c9
        dE[2] = 1.2 * c4 - 1.1 * c0
        for i in range(3):
            dE[i] = dE[i]*kappa*sqrt1_2

    else:
        s0 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        s4 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)
        sk = sin(kappa * x[2])
        ck = cos(kappa * x[2])

        dE[0] = 1.2 * s4 * sk + 1.3 * sqrt1_2 * c9 * ck
        dE[1] = -1.1 * s0 * sk - 1.3 * sqrt1_2 * c9 * ck
        dE[2] = -sqrt1_2 * (1.1 * c0 - 1.2 * c4) * ck
        for i in range(3):
            dE[i] = dE[i]*kappa


def f_exact(x, f):
    if dim == 1:
        s0 = sin(kappa * x[0] + 0.0 * pi)
        s4 = sin(kappa * x[0] + 0.4 * pi)
        s9 = sin(kappa * x[0] + 0.9 * pi)

        f[0] = 2.2 * s0 + 1.2 * sqrt1_2 * s4
        f[1] = (1.2 * (2.0 + kappa * kappa) * s4 +
                sqrt1_2 * (1.1 * s0 + 1.3 * s9))
        f[2] = 1.3 * (2.0 + kappa * kappa) * s9 + 1.2 * sqrt1_2 * s4

    elif dim == 2:
        s0 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        s4 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        s9 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)

        f[0] = (0.55 * (4.0 + kappa * kappa) * s0 +
                0.6 * (sqrt2 - kappa * kappa) * s4)
        f[1] = (0.55 * (sqrt2 - kappa * kappa) * s0 +
                0.6 * (4.0 + kappa * kappa) * s4 +
                0.65 * sqrt2 * s9)
        f[2] = 0.6 * sqrt2 * s4 + 1.3 * (2.0 + kappa * kappa) * s9

    else:
        s0 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        s4 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        s9 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)
        c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)
        sk = sin(kappa * x[2])
        ck = cos(kappa * x[2])

        f[0] = (0.55 * (4.0 + 3.0 * kappa * kappa) * s0 * ck +
                0.6 * (sqrt2 - kappa * kappa) * s4 * ck -
                0.65 * sqrt2 * kappa * kappa * c9 * sk)

        f[1] = (0.55 * (sqrt2 - kappa * kappa) * s0 * ck +
                0.6 * (4.0 + 3.0 * kappa * kappa) * s4 * ck +
                0.65 * sqrt2 * s9 * ck -
                0.65 * sqrt2 * kappa * kappa * c9 * sk)

        f[2] = (0.6 * sqrt2 * s4 * ck -
                sqrt2 * kappa * kappa * (0.55 * c0 + 0.6 * c4) * sk
                + 1.3 * (2.0 + kappa * kappa) * s9 * ck)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex31 (Definite Maxwell Problem)')
    parser.add_argument('-m', '--mesh',
                        default="inline-quad.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-r', '--refine',
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=1.0,
                        help="Set the frequency for the exact")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    try:
        from numba import jit
        HAS_NUMBA = True
    except ImportError:
        assert False, "This example requires numba to run"
    parser.add_argument("-n", "--numba",
                        default=int(HAS_NUMBA),
                        type=int,
                        help="Use Number compiled coefficient")

    args = parser.parse_args()
    args.numba = bool(args.numba)
    parser.print_options(args)

    order = args.order

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    freq = args.frequency
    numba = args.numba

    run(freq=freq,
        order=order,
        refine=args.refine,
        meshfile=meshfile,
        visualization=visualization,
        numba=numba)
