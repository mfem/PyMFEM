'''
   MFEM example 31

   See c++ version in the MFEM library for more detail

   Sample runs: mpirun -np 4 python ex31p.py -m ../data/hexagon.mesh -o 2
                mpirun -np 4 python ex31p.py -m ../data/star.mesh 
                mpirun -np 4 python ex31p.py -m ../data/square-disc.mesh -o 2
                mpirun -np 4 python ex31p.py -m ../data/fichera.mesh -o 3 -rs 1 -rp 0
                mpirun -np 4 python ex31p.py -m ../data/square-disc-nurbs.mesh -o 3
                mpirun -np 4 python ex31p.py -m ../data/amr-quad.mesh -o 2 -rs 1
                mpirun -np 4 python ex31p.py -m ../data/amr-hex.mesh -rs 1
'''
from mpi4py import MPI
import mfem.par as mfem
from mfem.par import intArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt

sqrt1_2 = 1/sqrt(2)
sqrt2 = sqrt(2)

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(order=1,
        rs=2,
        rp=1,
        freq=1,
        meshfile='',
        visualization=False,
        numba=False):

    # 3. Read the (serial) mesh from the given mesh file on all processors.  We
    #    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    #    and volume meshes with the same code.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    # 4. Refine the serial mesh on all processors to increase the resolution. In
    #    this example we do 'ref_levels' of uniform refinement (2 by default, or
    #    specified on the command line with -rs).
    for x in range(rs):
        mesh.UniformRefinement()

    # 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution (1 time by
    #    default, or specified on the command line with -rp). Once the parallel
    #    mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh
    for lev in range(rp):
        pmesh.UniformRefinement()

    # 6. Define a parallel finite element space on the parallel mesh. Here we
    #    use the Nedelec finite elements of the specified order.
    if dim == 1:
        fec = mfem.ND_R1D_FECollection(order, dim)
    elif dim == 2:
        fec = mfem.ND_R2D_FECollection(order, dim)
    else:
        fec = mfem.ND_FECollection(order, dim)

    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    size = fespace.GlobalTrueVSize()
    if myid == 0:
        print("Number of H(curl) unknowns: " + str(size))

    # 7. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, the boundary conditions are defined
    #    by marking all the boundary attributes from the mesh as essential
    #    (Dirichlet) and converting them to a list of true dofs.
    ess_tdof_list = intArray()
    if pmesh.bdr_attributes.Size():
        ess_bdr = intArray([1]*pmesh.bdr_attributes.Max())
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 7.b Define coefficent to use. Here, we define outside this function.
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

    # 8. Set up the parallel linear form b(.) which corresponds to the
    #    right-hand side of the FEM linear system, which in this case is
    #    (f,phi_i) where f is given by the function f_exact and phi_i are the
    #    basis functions in the finite element fespace.
    b = mfem.ParLinearForm(fespace)
    b.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(f))
    b.Assemble()

    # 9. Define the solution vector x as a parallel finite element grid function
    #    corresponding to fespace. Initialize x by projecting the exact
    #    solution. Note that only values from the boundary edges will be used
    #    when eliminating the non-homogeneous boundary condition to modify the
    #    r.h.s. vector b.
    sol = mfem.ParGridFunction(fespace)
    sol.ProjectCoefficient(E)

    # 10. Set up the parallel bilinear form corresponding to the EM diffusion
    #     operator curl muinv curl + sigma I, by adding the curl-curl and the
    #     mass domain integrators.
    mat = array([[2.0, sqrt1_2, 0, ],
                 [sqrt1_2, 2.0, sqrt1_2],
                 [0.0, sqrt1_2, 2.0, ], ])
    sigmaMat = mfem.DenseMatrix(mat)

    muinv = mfem.ConstantCoefficient(1.0)
    sigma = mfem. MatrixConstantCoefficient(sigmaMat)
    a = mfem.ParBilinearForm(fespace)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))

    # 11. Assemble the parallel bilinear form and the corresponding linear
    #     system, applying any necessary transformations such as: parallel
    #     assembly, eliminating boundary conditions, applying conforming
    # b    constraints for non-conforming AMR, etc.
    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, sol, b, A, X, B)

    # 12. Solve the system AX=B using PCG with the AMS preconditioner from hypre
    AM = A.AsHypreParMatrix()
    if myid == 0:
        print("Size of linear system: " + str(AM.GetGlobalNumRows()))

    ams = mfem.HypreAMS(AM, fespace)

    pcg = mfem.HyprePCG(AM)
    pcg.SetTol(1e-12)
    pcg.SetMaxIter(1000)
    pcg.SetPrintLevel(2)
    pcg.SetPreconditioner(ams)
    pcg.Mult(B, X)

    # 13. Recover the parallel grid function corresponding to X. This is the
    #     local finite element solution on each processor.
    a.RecoverFEMSolution(X, b, sol)

    # 14. Compute and print the CurlE norm of the error.
    error = sol.ComputeHCurlError(E, CurlE)
    if myid == 0:
        print("|| E_h - E ||_{Hcurl} = " + "{:g}".format(error)+"\n")

    # 15. Save the refined mesh and the solution in parallel. This output can
    #     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    smyid = '{:0>6d}'.format(myid)
    pmesh.Print('mesh.'+smyid, 8)
    sol.Save('sol.'+smyid, 8)

    # 15. Send the solution by socket to a GLVis server.
    if visualization:

        solCoef = mfem.VectorGridFunctionCoefficient(sol)
        dsolCoef = mfem.CurlGridFunctionCoefficient(sol)

        if dim == 1:
            x_sock = make_socketstrema()
            y_sock = make_socketstrema()
            z_sock = make_socketstrema()
            dy_sock = make_socketstrema()
            dz_sock = make_socketstrema()

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

            fes_h1 = mfem.ParFiniteElementSpace(pmesh, fec_h1)
            fes_l2 = mfem.ParFiniteElementSpace(pmesh, fec_l2)

            xComp = mfem.ParGridFunction(fes_l2)
            yComp = mfem.ParGridFunction(fes_h1)
            zComp = mfem.ParGridFunction(fes_h1)

            dyComp = mfem.ParGridFunction(fes_l2)
            dzComp = mfem.ParGridFunction(fes_l2)

            xCoef = mfem.InnerProductCoefficient(xVecCoef, solCoef)
            yCoef = mfem.InnerProductCoefficient(yVecCoef, solCoef)
            zCoef = mfem.InnerProductCoefficient(zVecCoef, solCoef)

            xComp.ProjectCoefficient(xCoef)
            yComp.ProjectCoefficient(yCoef)
            zComp.ProjectCoefficient(zCoef)

            x_sock << "solution\n" << pmesh << xComp
            x_sock.flush()
            x_sock << "window_title 'X component'"
            x_sock.flush()
            y_sock << "solution\n" << pmesh << yComp
            y_sock.flush()
            y_sock << "window_geometry 403 0 400 350 " << "window_title 'Y component'"
            y_sock.flush()
            z_sock << "solution\n" << pmesh << zComp
            z_sock.flush()
            z_sock << "window_geometry 806 0 400 350 " << "window_title 'Z component'"
            z_sock.flush()

            dyCoef = mfem.InnerProductCoefficient(yVecCoef, dsolCoef)
            dzCoef = mfem.InnerProductCoefficient(zVecCoef, dsolCoef)
            dyComp.ProjectCoefficient(dyCoef)
            dzComp.ProjectCoefficient(dzCoef)

            dy_sock << "solution\n" << pmesh << dyComp
            dy_sock.flush()
            dy_sock << "window_geometry 403 375 400 350 " << "window_title 'Y component of Curl'"
            dy_sock.flush()

            dy_sock << "solution\n" << pmesh << dzComp
            dy_sock.flush()
            dy_sock << "window_geometry 403 375 400 350 " << "window_title 'Z component of Curl'"
            dy_sock.flush()

        elif dim == 2:
            xy_sock = make_socketstrema()
            z_sock = make_socketstrema()
            dxy_sock = make_socketstrema()
            dz_sock = make_socketstrema()

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

            fes_h1 = mfem.ParFiniteElementSpace(pmesh, fec_h1)
            fes_nd = mfem.ParFiniteElementSpace(pmesh, fec_nd)
            fes_rt = mfem.ParFiniteElementSpace(pmesh, fec_rt)
            fes_l2 = mfem.ParFiniteElementSpace(pmesh, fec_l2)

            xyComp = mfem.ParGridFunction(fes_nd)
            zComp = mfem.ParGridFunction(fes_h1)

            dxyComp = mfem.ParGridFunction(fes_rt)
            dzComp = mfem.ParGridFunction(fes_l2)

            xyComp.ProjectCoefficient(xyCoef)
            zComp.ProjectCoefficient(zCoef)

            xy_sock.precision(8)
            xy_sock << "solution\n" << pmesh << xyComp
            xy_sock << "window_title 'XY components'\n"
            xy_sock.flush()
            z_sock << "solution\n" << pmesh << zComp
            z_sock << "window_geometry 403 0 400 350 " << "window_title 'Z component'"
            z_sock.flush()

            dxyCoef = mfem.MatrixVectorProductCoefficient(xyMatCoef, dsolCoef)
            dzCoef = mfem.InnerProductCoefficient(zVecCoef, dsolCoef)

            dxyComp.ProjectCoefficient(dxyCoef)
            dzComp.ProjectCoefficient(dzCoef)

            dxy_sock << "solution\n" << pmesh << dxyComp
            dxy_sock << "window_geometry 0 375 400 350 " << "window_title 'XY components of Curl'"
            dxy_sock.flush()
            dz_sock << "solution\n" << pmesh << dzComp
            dz_sock << "window_geometry 403 375 400 350 " << "window_title 'Z component of Curl'"
            dz_sock.flush()

        else:
            sol_sock = make_socketstrema()
            dsol_sock = make_socketstrema()

            fec_rt = mfem.RT_FECollection(order-1, dim)
            fes_rt = mfem.ParFiniteElementSpace(pmesh, fec_rt)

            dsol = mfem.ParGridFunction(fes_rt)

            dsol.ProjectCoefficient(dsolCoef)

            sol_sock.precision(8)
            dsol_sock.precision(8)
            sol_sock << "solution\n" << pmesh << sol << "window_title 'Solution'"
            sol_sock.flush()
            dsol_sock << "solution\n" << pmesh << dsol
            dsol_sock.flush()
            dsol_sock << "window_geometry 0 375 400 350 " << "window_title 'Curl of solution'"
            dsol_sock.flush()


def make_socketstrema():
    sock = mfem.socketstream("localhost", 19916)
    sock.precision(8)
    sock << "parallel " << num_procs << " " << myid << "\n"
    return sock


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
    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly in serial")
    parser.add_argument('-rp', '--refine-parallel',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly in paralle.")
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
    if myid == 0:
        parser.print_options(args)

    order = args.order

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    freq = args.frequency
    numba = args.numba

    run(freq=freq,
        order=order,
        rs=args.refine_serial,
        rp=args.refine_parallel,
        meshfile=meshfile,
        visualization=visualization,
        numba=numba)
