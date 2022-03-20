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
import mfem.par as mfem
from mfem.par import intArray, doubleArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt

sqrt1_2 = 1/sqrt(2)
sqrt2 = sqrt(2)

from mpi4py import MPI
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
        print("Number of H(Div) unknowns: " + str(size_nd))

    # 7. Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
    #    element space. The first corresponds to the curl curl, while the second
    #    is a simple mass matrix needed on the right hand side of the
    #    generalized eigenvalue problem below. The boundary conditions are
    #    implemented by marking all the boundary attributes from the mesh as
    #    essential. The corresponding degrees of freedom are eliminated with
    #    special values on the diagonal to shift the Dirichlet eigenvalues out
    #    of the computational range. After serial and parallel assembly we
    #    extract the corresponding parallel matrices A and M.        

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
         shift = 1.0;
         if myid == 0:
            print("Computing eigenvalues shifted by " + str(1.0))

    a.Assemble();
    a.EliminateEssentialBCDiag(ess_bdr, 1.0);
    a.Finalize()

    m = mfem.ParBilinearForm(fespace_nd);
    m.AddDomainIntegrator(mass)
    m.Assemble();

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
    curl.AddDomainInterpolator(mfem.CurlInterpolator());
    curl.Assemble();
    curl.Finalize();

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

    # 13. Send the solution by socket to a GLVis server.    
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

            xVec = mdfem.Vector([1, 0, 0.])
            yVec = mdfem.Vector([0, 1, 0.])
            zVec = mdfem.Vector([0, 0, 1.])

            xVecCof = mfem.VectorConstantCoefficient(xVec)
            yVecCof = mfem.VectorConstantCoefficient(yVec)
            zVecCof = mfem.VectorConstantCoefficient(zVec)

            fec_h1 = mfem.H1_FECollection(order, dim)
            fec_l2 = mfem.L2_FECollection(order-1, dim)

            fes_h1 = mfem.FiniteElementSpace(mesh, fec_h1)
            fes_l2 = mfem.FiniteElementSpace(mesh, fec_l2)

            xComp = mfem.GridFunction(fes_l2)
            yComp = mfem.GridFunction(fes_h1)
            zComp = mfem.GridFunction(fes_h1)

            xCoef = mfem.InnerProductCoefficient(xVecCoef, solCoef)
            yCoef = mfem.InnerProductCoefficient(yVecCoef, solCoef)
            zCoef = mfem.InnerProductCoefficient(zVecCoef, solCoef)

            xComp.ProjectCoefficient(xCoef)
            yComp.ProjectCoefficient(yCoef)
            zComp.ProjectCoefficient(zCoef)

            x_sock << "solution\n" << mesh << xComp
            x_sock.fluxh()
            x_sock << "window_title 'X component'"
            x_sock.fluxh()
            y_sock << "solution\n" << mesh << yComp
            y_sock.fluxh()
            y_sock << "window_geometry 403 0 400 350 " << "window_title 'Y component'"
            y_sock.fluxh()
            z_sock << "solution\n" << mesh << zComp
            z_sock.fluxh()
            z_sock << "window_geometry 806 0 400 350 " << "window_title 'Z component'"
            z_sock.fluxh()

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

            zVec = array([0, 0, 1])
            zVecCoef = mfem.VectorConstantCoefficient(zVec)

            xyCoeff = mfem.MatrixVectorProductCoefficient(xyMatCoef, solCoef)
            zCeof = mfem.InnerProductCoefficient(zVecCoef, solCoef)

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
            zComp = mfem.GridFunction(fes_l2)

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
        
        for i in range(3):  E[i] = E[i]*cos(kappa * x[2])

def CurlE_exact(x, dE):
    if dim == 1:
        c4 = cos(kappa * x[0] + 0.4 * pi)
        c9 = cos(kappa * x[0] + 0.9 * pi)

        dE[0] = 0.0
        dE[1] = -1.3 * c9
        dE[2] = 1.2 * c4
        for i in range(3):  dE[i] = dE[i]*kappa
        
    elif dim == 2:
        c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
        c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
        c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)

        dE[0] = 1.3 * c9
        dE[1] = -1.3 * c9
        dE[2] = 1.2 * c4 - 1.1 * c0
        for i in range(3):  dE[i] = dE[i]*kappa*sqrt1_2        

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
        for i in range(3):  dE[i] = dE[i]*kappa               

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
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization

    run(nev=args.num_eigs,
        order=order,
        rs=args.refine_serial,
        rp=args.refine_parallel,
        meshfile=meshfile,
        visualization=visualization)

        
