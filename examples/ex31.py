'''
   MFEM example 31

   See c++ version in the MFEM library for more detail 
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

    kappa = pi*freq

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
        params = {"dim":dim}
        E = mfem.jit.vector(E_exact, sdim=3, params=params)
        CurlE = mfem.jit.vector(curlE_exact, sdim=3, params=params)
        f = mfem.jit.vector(f_exact, sdim=3, params=params)        
    else:
        pass # ToDo provide regular example.

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
    sigmaMat = mfem.DenseMatrix(3);
    sigmaMat[0,0] = 2.0
    sigmaMat[1,1] = 2.0
    sigmaMat[2,2] = 2.0
    sigmaMat[0,2] = 0.0
    sigmaMat[2,0] = 0.0
    sigmaMat[0,1] = sqrt1_2
    sigmaMat[1,0] = sqrt1_2
    sigmaMat[1,2] = sqrt1_2
    sigmaMat[2,1] = sqrt1_2

    muinv = mfem.ConstantCoefficient(1.0)
    sigma = mfem. MatrixConstantCoefficient(sigmaMat)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))
    
    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
    

    # 10. Solve the system A X = B.
    M = mfem.GSSmoother(A.Ptr())
    mfem.PCG(A.Ptr(), M, B, X, 1, 500, 1e-12, 0.0)

    # 12. Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, x)

    # 13. Compute and print the CurlE norm of the error.
    print("|| E_h - E ||_{Hcurl} = " + str(x.ComputeCurlError(CurlE_exact))+"\n")

    # 14. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g sol.gf"    
    mesh.Print('refined.mesh', 8)
    x.Save('sol.gf', 8)


// 15. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;

      VectorGridFunctionCoefficient solCoef(&sol);
      CurlGridFunctionCoefficient dsolCoef(&sol);

      if (dim ==1)
      {
         socketstream x_sock(vishost, visport);
         socketstream y_sock(vishost, visport);
         socketstream z_sock(vishost, visport);
         socketstream dy_sock(vishost, visport);
         socketstream dz_sock(vishost, visport);
         x_sock.precision(8);
         y_sock.precision(8);
         z_sock.precision(8);
         dy_sock.precision(8);
         dz_sock.precision(8);

         Vector xVec(3); xVec = 0.0; xVec(0) = 1;
         Vector yVec(3); yVec = 0.0; yVec(1) = 1;
         Vector zVec(3); zVec = 0.0; zVec(2) = 1;
         VectorConstantCoefficient xVecCoef(xVec);
         VectorConstantCoefficient yVecCoef(yVec);
         VectorConstantCoefficient zVecCoef(zVec);

         H1_FECollection fec_h1(order, dim);
         L2_FECollection fec_l2(order-1, dim);

         FiniteElementSpace fes_h1(&mesh, &fec_h1);
         FiniteElementSpace fes_l2(&mesh, &fec_l2);

         GridFunction xComp(&fes_l2);
         GridFunction yComp(&fes_h1);
         GridFunction zComp(&fes_h1);

         GridFunction dyComp(&fes_l2);
         GridFunction dzComp(&fes_l2);

         InnerProductCoefficient xCoef(xVecCoef, solCoef);
         InnerProductCoefficient yCoef(yVecCoef, solCoef);
         InnerProductCoefficient zCoef(zVecCoef, solCoef);

         xComp.ProjectCoefficient(xCoef);
         yComp.ProjectCoefficient(yCoef);
         zComp.ProjectCoefficient(zCoef);

         x_sock << "solution\n" << mesh << xComp << flush
                << "window_title 'X component'" << endl;
         y_sock << "solution\n" << mesh << yComp << flush
                << "window_geometry 403 0 400 350 "
                << "window_title 'Y component'" << endl;
         z_sock << "solution\n" << mesh << zComp << flush
                << "window_geometry 806 0 400 350 "
                << "window_title 'Z component'" << endl;

         InnerProductCoefficient dyCoef(yVecCoef, dsolCoef);
         InnerProductCoefficient dzCoef(zVecCoef, dsolCoef);

         dyComp.ProjectCoefficient(dyCoef);
         dzComp.ProjectCoefficient(dzCoef);

         dy_sock << "solution\n" << mesh << dyComp << flush
                 << "window_geometry 403 375 400 350 "
                 << "window_title 'Y component of Curl'" << endl;
         dz_sock << "solution\n" << mesh << dzComp << flush
                 << "window_geometry 806 375 400 350 "
                 << "window_title 'Z component of Curl'" << endl;
      }
      else if (dim == 2)
      {
         socketstream xy_sock(vishost, visport);
         socketstream z_sock(vishost, visport);
         socketstream dxy_sock(vishost, visport);
         socketstream dz_sock(vishost, visport);

         DenseMatrix xyMat(2,3); xyMat = 0.0;
         xyMat(0,0) = 1.0; xyMat(1,1) = 1.0;
         MatrixConstantCoefficient xyMatCoef(xyMat);
         Vector zVec(3); zVec = 0.0; zVec(2) = 1;
         VectorConstantCoefficient zVecCoef(zVec);

         MatrixVectorProductCoefficient xyCoef(xyMatCoef, solCoef);
         InnerProductCoefficient zCoef(zVecCoef, solCoef);

         H1_FECollection fec_h1(order, dim);
         ND_FECollection fec_nd(order, dim);
         RT_FECollection fec_rt(order-1, dim);
         L2_FECollection fec_l2(order-1, dim);

         FiniteElementSpace fes_h1(&mesh, &fec_h1);
         FiniteElementSpace fes_nd(&mesh, &fec_nd);
         FiniteElementSpace fes_rt(&mesh, &fec_rt);
         FiniteElementSpace fes_l2(&mesh, &fec_l2);

         GridFunction xyComp(&fes_nd);
         GridFunction zComp(&fes_h1);

         GridFunction dxyComp(&fes_rt);
         GridFunction dzComp(&fes_l2);

         xyComp.ProjectCoefficient(xyCoef);
         zComp.ProjectCoefficient(zCoef);

         xy_sock.precision(8);
         xy_sock << "solution\n" << mesh << xyComp
                 << "window_title 'XY components'\n" << flush;
         z_sock << "solution\n" << mesh << zComp << flush
                << "window_geometry 403 0 400 350 "
                << "window_title 'Z component'" << endl;

         MatrixVectorProductCoefficient dxyCoef(xyMatCoef, dsolCoef);
         InnerProductCoefficient dzCoef(zVecCoef, dsolCoef);

         dxyComp.ProjectCoefficient(dxyCoef);
         dzComp.ProjectCoefficient(dzCoef);

         dxy_sock << "solution\n" << mesh << dxyComp << flush
                  << "window_geometry 0 375 400 350 "
                  << "window_title 'XY components of Curl'" << endl;
         dz_sock << "solution\n" << mesh << dzComp << flush
                 << "window_geometry 403 375 400 350 "
                 << "window_title 'Z component of Curl'" << endl;
      }
      else
      {
         socketstream sol_sock(vishost, visport);
         socketstream dsol_sock(vishost, visport);

         RT_FECollection fec_rt(order-1, dim);

         FiniteElementSpace fes_rt(&mesh, &fec_rt);

         GridFunction dsol(&fes_rt);

         dsol.ProjectCoefficient(dsolCoef);

         sol_sock.precision(8);
         sol_sock << "solution\n" << mesh << sol
                  << "window_title 'Solution'" << flush << endl;
         dsol_sock << "solution\n" << mesh << dsol << flush
                   << "window_geometry 0 375 400 350 "
                   << "window_title 'Curl of solution'" << endl;
      }
   }
    

if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex3 (Definite Maxwell Problem)')
    parser.add_argument('-m', '--mesh',
                        default="beam-tet.mesh",
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
    refine = args.refine

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    freq = args.frequency
    numba = args.numba

    run(freq=freq,
        order=order,
        refine=refine,
        meshfile=meshfile,
        visualization=visualization,
        numba=numba)

def E_exact(x, out):
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
       E *= cos(kappa * x[2])


def CurlE_exact(x, dE):
    if dim == 1:
      c4 = cos(kappa * x[0] + 0.4 * pi)
      c9 = cos(kappa * x[0] + 0.9 * pi)

      dE[0] =  0.0
      dE[1] = -1.3 * c9
      dE[2] =  1.2 * c4
      dE *= kappa;
    elif dim == 2:
      c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
      c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
      c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)

      dE[0] =  1.3 * c9
      dE[1] = -1.3 * c9
      dE[2] =  1.2 * c4 - 1.1 * c0
      dE *= kappa * sqrt1_2;
    else:
       s0 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
       c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi)
       s4 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
       c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi)
       c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi)
       sk = sin(kappa * x[2])
       ck = cos(kappa * x[2])
 
       dE[0] =  1.2 * s4 * sk + 1.3 * sqrt1_2 * c9 * ck
       dE[1] = -1.1 * s0 * sk - 1.3 * sqrt1_2 * c9 * ck
       dE[2] = -sqrt1_2 * (1.1 * c0 - 1.2 * c4) * ck
       dE *= kappa;

def f_exact(x, f):
   if dim == 1:
      s0 = sin(kappa * x[0] + 0.0 * pi)
      s4 = sin(kappa * x[0] + 0.4 * pi)
      s9 = sin(kappa * x[0] + 0.9 * pi)

      f[0] = 2.2 * s0 + 1.2 * sqrt1_2 * s4;
      f[1] = 1.2 * (2.0 + kappa * kappa) * s4 +
             sqrt1_2 * (1.1 * s0 + 1.3 * s9);
      f[2] = 1.3 * (2.0 + kappa * kappa) * s9 + 1.2 * sqrt1_2 * s4;
      
   elif dim == 2:
      double s0 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi);
      double s4 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi);
      double s9 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi);

      f[0] = 0.55 * (4.0 + kappa * kappa) * s0 +
             0.6 * (sqrt2 - kappa * kappa) * s4;
      f[1] = 0.55 * (sqrt2 - kappa * kappa) * s0 +
             0.6 * (4.0 + kappa * kappa) * s4 +
             0.65 * sqrt2 * s9;
      f[2] = 0.6 * sqrt2 * s4 + 1.3 * (2.0 + kappa * kappa) * s9;

   else:

      double s0 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi);
      double c0 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.0 * pi);
      double s4 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi);
      double c4 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.4 * pi);
      double s9 = sin(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi);
      double c9 = cos(kappa * sqrt1_2 * (x[0] + x[1]) + 0.9 * pi);
      double sk = sin(kappa * x[2]);
      double ck = cos(kappa * x[2]);

      f[0] = 0.55 * (4.0 + 3.0 * kappa * kappa) * s0 * ck +
             0.6 * (sqrt2 - kappa * kappa) * s4 * ck -
             0.65 * sqrt2 * kappa * kappa * c9 * sk;

      f[1] = 0.55 * (sqrt2 - kappa * kappa) * s0 * ck +
             0.6 * (4.0 + 3.0 * kappa * kappa) * s4 * ck +
             0.65 * sqrt2 * s9 * ck -
             0.65 * sqrt2 * kappa * kappa * c9 * sk;

      f[2] = 0.6 * sqrt2 * s4 * ck -
             sqrt2 * kappa * kappa * (0.55 * c0 + 0.6 * c4) * sk
             + 1.3 * (2.0 + kappa * kappa) * s9 * ck;

