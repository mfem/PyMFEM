'''
   MFEM example 32

   See c++ version in the MFEM library for more detail

   Sample runs:  mpirun -np 4 python ex32p.py -m ../data/hexagon.mesh -o 2
                 mpirun -np 4 python ex32p.py -m ../data/star.mesh
                 mpirun -np 4 python ex32p.py -m ../data/square-disc.mesh -o 2 -n 4 -rs 1
                 mpirun -np 4 python ex32p.py -m ../data/square-disc-nurbs.mesh -rs 3 -o 3
                 mpirun -np 4 python ex32p.py -m ../data/amr-quad.mesh -o 2 -rs 1
                 mpirun -np 4 python ex32p.py -m ../data/amr-hex.mesh -rs 1
                 mpirun -np 4 python ex32p.py -m ../data/fichera.mesh -rs 1

'''
from mpi4py import MPI
import mfem.par as mfem
from mfem.par import intArray, doubleArray
import os
import sys
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt

sqrt1_2 = 1/sqrt(2)
sqrt2 = sqrt(2)

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def run(nev=5,
        order=1,
        rs=2,
        rp=1,
        meshfile='',
        visualization=False):

    # 3. Read the (serial) mesh from the given mesh file on all processors. We
    #    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    #    and volume meshes with the same code.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    bbMin, bbMax = mesh.GetBoundingBox()

    # 4. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement (2 by default, or specified on
    #     the command line with -rs)
    for lev in range(rs):
        mesh.UniformRefinement()

    # 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution (1 time by
    #    default, or specified on the command line with -rp). Once the parallel

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh
    for lev in range(rp):
        pmesh.UniformRefinement()

    # 6. Define a parallel finite element space on the parallel mesh. Here we
    #    use the Nedelec finite elements of the specified order.
    if dim == 1:
        fec_nd = mfem.ND_R1D_FECollection(order, dim)
        fec_rt = mfem.RT_R1D_FECollection(order-1, dim)
    elif dim == 2:
        fec_nd = mfem.ND_R2D_FECollection(order, dim)
        fec_rt = mfem.RT_R2D_FECollection(order-1, dim)
    else:
        fec_nd = mfem.ND_FECollection(order, dim)
        fec_rt = mfem.RT_FECollection(order-1, dim)

    fespace_nd = mfem.ParFiniteElementSpace(pmesh, fec_nd)
    fespace_rt = mfem.ParFiniteElementSpace(pmesh, fec_rt)

    size_nd = fespace_nd.GlobalTrueVSize()
    size_rt = fespace_rt.GlobalTrueVSize()

    if myid == 0:
        print("Number of H(Curl) unknowns: " + str(size_nd))
        print("Number of H(Div) unknowns: " + str(size_rt))

    # 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
    #    element space. The first corresponds to the curl curl, while the second
    #    is a simple mass matrix needed on the right hand side of the
    #    generalized eigenvalue problem below. The boundary conditions are
    #    implemented by marking all the boundary attributes from the mesh as
    #    essential. The corresponding degrees of freedom are eliminated with
    #    special values on the diagonal to shift the Dirichlet eigenvalues out
    #    of the computational range. After serial and parallel assembly we
    #    extract the corresponding parallel matrices A and M.
    shift = 0.0
    mat = array([[2.0, sqrt1_2, 0, ],
                 [sqrt1_2, 2.0, sqrt1_2],
                 [0.0, sqrt1_2, 2.0, ], ])
    epsilon_mat = mfem.DenseMatrix(mat)
    epsilon = mfem. MatrixConstantCoefficient(epsilon_mat)
    one = mfem.ConstantCoefficient(1.0)

    if pmesh.bdr_attributes.Size():
        ess_bdr = intArray([1]*pmesh.bdr_attributes.Max())

    a = mfem.ParBilinearForm(fespace_nd)
    curlcurl = mfem.CurlCurlIntegrator(one)
    a.AddDomainIntegrator(curlcurl)

    mass = mfem.VectorFEMassIntegrator(epsilon)

    if pmesh.bdr_attributes.Size() == 0 or dim == 1:
        #  Add a mass term if the mesh has no boundary, e.g. periodic mesh or
        #  closed surface.
        a.AddDomainIntegrator(mass)
        shift = 1.0
        if myid == 0:
            print("Computing eigenvalues shifted by " + str(1.0))

    a.Assemble()
    a.EliminateEssentialBCDiag(ess_bdr, 1.0)
    a.Finalize()

    m = mfem.ParBilinearForm(fespace_nd)
    m.AddDomainIntegrator(mass)
    m.Assemble()

    # shift the eigenvalue corresponding to eliminated dofs to a large value
    m.EliminateEssentialBCDiag(ess_bdr, 2.3e-308)
    m.Finalize()

    A = a.ParallelAssemble()
    M = m.ParallelAssemble()

    # 8. Define and configure the AME eigensolver and the AMS preconditioner for
    #    A to be used within the solver. Set the matrices which define the
    #    generalized eigenproblem A x = lambda M x.
    ams = mfem.HypreAMS(A, fespace_nd)
    ams.SetPrintLevel(0)
    ams.SetSingularProblem()

    ame = mfem.HypreAME(MPI.COMM_WORLD)
    ame.SetNumModes(nev)
    ame.SetPreconditioner(ams)
    ame.SetMaxIter(100)
    ame.SetTol(1e-8)
    ame.SetPrintLevel(1)
    ame.SetMassMatrix(M)
    ame.SetOperator(A)

    # 9. Compute the eigenmodes and extract the array of eigenvalues. Define
    #    parallel grid functions to represent each of the eigenmodes returned by
    #    the solver and their derivatives.
    eigenvalues = doubleArray()
    ame.Solve()
    ame.GetEigenvalues(eigenvalues)
    x = mfem.ParGridFunction(fespace_nd)
    dx = mfem.ParGridFunction(fespace_rt)

    curl = mfem.ParDiscreteLinearOperator(fespace_nd, fespace_rt)
    curl.AddDomainInterpolator(mfem.CurlInterpolator())
    curl.Assemble()
    curl.Finalize()

    # 10. Save the refined mesh and the modes in parallel. This output can be
    #     viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
    smyid = '{:0>6d}'.format(myid)
    mesh_name = "mesh."+smyid
    pmesh.Print(mesh_name, 8)

    for i in range(nev):
        x.Assign(ame.GetEigenvector(i))
        curl.Mult(x, dx)

        mode_name = "mode_"+str(i).zfill(2)+"."+smyid
        mode_deriv_name = "mode_deriv_"+str(i).zfill(2)+"."+smyid

        x.Save(mode_name, 8)
        dx.Save(mode_deriv_name, 8)

    # 11 Visualize data using glvis
    if visualization:
        # 11. (a) functions to format Glvis window
        def make_cmd(label, px, py, max_r, i, title="Eigenmode"):
            cmd = (" window_title '" + title + " " + str(i+1) + '/' + str(nev) +
                   " " + label + ", Lambda = " + "{:g}".format(eigenvalues[i] - shift) +
                   "'" + " valuerange -" + str(max_r) + ' ' + str(max_r))
            cmd = (cmd + " keys aa window_geometry " +
                   str(px) + " " + str(py) + " 400 350")
            return cmd
        # 11. (b) fucntion to send data to Glvis

        def send_data_to_glvis(sock, data, cmd):
            sock << "parallel " << num_procs << " " << myid << "\n"
            sock << "solution\n" << pmesh << data
            sock.flush()
            sock << cmd
            sock.endline()
            MPI.COMM_WORLD.Barrier()

        def ask_exit():
            if (myid == 0):
                from builtins import input
                c = input("press (q)uit or (c)ontinue --> ")
            else:
                c = None
            c = MPI.COMM_WORLD.bcast(c, root=0)
            if (c != 'c'):
                sys.exit()

        if dim == 1:
            mode_x_sock = make_socketstrema()
            mode_y_sock = make_socketstrema()
            mode_z_sock = make_socketstrema()
            mode_dy_sock = make_socketstrema()
            mode_dz_sock = make_socketstrema()

            xVec = mdfem.Vector([1, 0, 0.])
            yVec = mdfem.Vector([0, 1, 0.])
            zVec = mdfem.Vector([0, 0, 1.])

            xVecCof = mfem.VectorConstantCoefficient(xVec)
            yVecCof = mfem.VectorConstantCoefficient(yVec)
            zVecCof = mfem.VectorConstantCoefficient(zVec)

            fec_h1 = mfem.H1_FECollection(order, dim)
            fec_l2 = mfem.L2_FECollection(order-1, dim)

            fes_h1 = mfem.ParFiniteElementSpace(pmesh, fec_h1)
            fes_l2 = mfem.ParFiniteElementSpace(pmesh, fec_l2)

            xComp = mfem.ParGridFunction(fes_l2)
            yComp = mfem.ParGridFunction(fes_h1)
            zComp = mfem.ParGridFunction(fes_h1)

            dyComp = mfem.ParGridFunction(fes_l2)
            dzComp = mfem.ParGridFunction(fes_l2)

            for i in range(nev):
                if (myid == 0):
                    print("Eigenmode " + str(i+1) + '/' + str(nev) +
                          ", Lambda = " + "{:g}".format(eigenvalues[i]))

                x.Assign(ame.GetEigenvector(i))
                curl.Mult(x, dx)

                modeCoeff = mfem.VectorGridFunctionCoefficient(x)
                xCoef = mfem.InnerProductCoefficient(xVecCoef, modeCoef)
                yCoef = mfem.InnerProductCoefficient(yVecCoef, modeCoef)
                zCoef = mfem.InnerProductCoefficient(zVecCoef, modeCoef)

                xComp.ProjectCoefficient(xCoef)
                yComp.ProjectCoefficient(yCoef)
                zComp.ProjectCoefficient(zCoef)

                max_x = GetScalarMax(xComp)
                max_y = GetScalarMax(yComp)
                max_z = GetScalarMax(zComp)
                max_r = np.max([max_x, max_y, max_z])

                x_cmd = make_cmd('X', 0, 0, max_r, i)
                y_cmd = make_cmd('Y', 403, 0, max_r, i)
                z_cmd = make_cmd('Z', 800, 0, max_r, i)

                send_data_to_glvis(mode_x_sock, x, x_cmd)
                send_data_to_glvis(mode_y_sock, y, y_cmd)
                send_data_to_glvis(mode_z_sock, z, z_cmd)

                dmodeCoeff = mfem.VectorGridFunctionCoefficient(dx)
                dyCoef = mfem.InnerProductCoefficient(yVecCoef, dmodeCoef)
                dzCoef = mfem.InnerProductCoefficient(zVecCoef, dmodeCoef)
                zCoef = mfem.InnerProductCoefficient(zVecCoef, modeCoef)

                dyComp.ProjectCoefficient(dyCoef)
                dzComp.ProjectCoefficient(dzCoef)

                min_d = max_r / bbMax[0] - bbMin[0]
                max_y = GetScalarMax(dyComp)
                max_z = GetScalarMax(dzComp)
                max_r = np.max([max_y, max_z, max_d])

                dy_cmd = make_cmd('Y', 403, 375, max_r, i,
                                  title="Curl Eigenmode")
                dz_cmd = make_cmd('Z', 800, 375, max_r, i,
                                  title="Curl Eigenmode")

                send_data_to_glvis(mode_dy_sock, dyComp, dy_cmd)
                send_data_to_glvis(mode_dz_sock, dzComp, dz_cmd)

                ask_exit()

        elif dim == 2:
            mode_xy_sock = make_socketstrema()
            mode_z_sock = make_socketstrema()
            mode_dxy_sock = make_socketstrema()
            mode_dz_sock = make_socketstrema()

            mat = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
            xyMat = mfem.DenseMatrix(mat)
            xyMatCoef = mfem.MatrixConstantCoefficient(xyMat)
            zVec = array([0, 0, 1])
            zVecCoef = mfem.VectorConstantCoefficient(zVec)

            fec_h1 = mfem.H1_FECollection(order, dim)
            fec_nd_xy = mfem.ND_FECollection(order, dim)
            fec_rt_xy = mfem.RT_FECollection(order-1, dim)
            fec_l2 = mfem.L2_FECollection(order-1, dim)

            fes_h1 = mfem.ParFiniteElementSpace(pmesh, fec_h1)
            fes_nd = mfem.ParFiniteElementSpace(pmesh, fec_nd_xy)
            fes_rt = mfem.ParFiniteElementSpace(pmesh, fec_rt_xy)
            fes_l2 = mfem.ParFiniteElementSpace(pmesh, fec_l2)

            xyComp = mfem.ParGridFunction(fes_nd)
            zComp = mfem.ParGridFunction(fes_h1)

            dxyComp = mfem.ParGridFunction(fes_rt)
            dzComp = mfem.ParGridFunction(fes_l2)

            for i in range(nev):
                if (myid == 0):
                    print("Eigenmode " + str(i+1) + '/' + str(nev) +
                          ", Lambda = " + "{:g}".format(eigenvalues[i]))

                x.Assign(ame.GetEigenvector(i))
                curl.Mult(x, dx)

                modeCoef = mfem.VectorGridFunctionCoefficient(x)

                xyCoef = mfem.MatrixVectorProductCoefficient(
                    xyMatCoef, modeCoef)
                zCoef = mfem.InnerProductCoefficient(zVecCoef, modeCoef)

                xyComp.ProjectCoefficient(xyCoef)
                zComp.ProjectCoefficient(zCoef)

                max_v = GetVectorMax(2, xyComp)
                max_s = GetScalarMax(zComp)
                max_r = np.max([max_v, max_s])

                xy_cmd = make_cmd('XY', 0, 0, max_r, i)
                z_cmd = make_cmd('Z', 403, 0, max_r, i)

                send_data_to_glvis(mode_xy_sock, xyComp, xy_cmd)
                send_data_to_glvis(mode_z_sock, zComp, z_cmd)

                dmodeCoef = mfem.VectorGridFunctionCoefficient(dx)
                dxyCoef = mfem.MatrixVectorProductCoefficient(
                    xyMatCoef, dmodeCoef)
                dzCoef = mfem.InnerProductCoefficient(zVecCoef, dmodeCoef)

                dxyComp.ProjectCoefficient(dxyCoef)
                dzComp.ProjectCoefficient(dzCoef)

                min_d = max_r / \
                    min([bbMax[0] - bbMin[0],  bbMax[1] - bbMin[1]])
                max_v = GetVectorMax(2, dxyComp)
                max_s = GetScalarMax(dzComp)
                max_r = np.max([max_v, max_s])

                dxy_cmd = make_cmd('XY', 0, 375, max_r, i,
                                   title="Curl Eigenmode")
                dz_cmd = make_cmd('Z', 403, 375, max_r, i,
                                  title="Curl Eigenmode")

                send_data_to_glvis(mode_dxy_sock, dxyComp, dxy_cmd)
                send_data_to_glvis(mode_dz_sock, dzComp, dz_cmd)

                ask_exit()

        else:
            mode_sock = make_socketstrema()
            dmode_sock = make_socketstrema()

            for i in range(nev):
                if (myid == 0):
                    print("Eigenmode " + str(i+1) + '/' + str(nev) +
                          ", Lambda = " + "{:g}".format(eigenvalues[i]))

                x.Assign(ame.GetEigenvector(i))
                curl.Mult(x, dx)

                max_v = GetVectorMax(3, x)
                cmd = make_cmd('', 0, 0, max_v, i)
                max_v = GetVectorMax(3, dx)
                d_cmd = make_cmd('', 0, 375, max_v, i)

                send_data_to_glvis(mode_sock, x, cmd)
                send_data_to_glvis(dmode_sock, dx, d_cmd)

                ask_exit()


def GetVectorMax(vdim, x):
    zeroVec = mfem.Vector(vdim)
    zeroVec.Assign(0.0)
    zero = mfem.VectorConstantCoefficient(zeroVec)
    nrm = x.ComputeMaxError(zero)
    return nrm


def GetScalarMax(x):
    zero = mfem.ConstantCoefficient(0.0)
    nrm = x.ComputeMaxError(zero)
    return nrm


def make_socketstrema():
    sock = mfem.socketstream("localhost", 19916)
    sock.precision(8)
    return sock


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex32p (Maxwell Eigenvalue Problem)')
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
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")
    parser.add_argument("-n", "--num-eigs",
                        action='store',
                        type=float,
                        default=5,
                        help="Number of eigen values to compute")
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = not args.no_visualization

    run(nev=args.num_eigs,
        order=order,
        rs=args.refine_serial,
        rp=args.refine_parallel,
        meshfile=meshfile,
        visualization=visualization)
