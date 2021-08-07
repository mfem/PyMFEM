'''
   MFEM example 29
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array


def GetMesh(mesh_type):
    if type == 3:
        mesh = mfem.Mesh(2, 12, 16, 8, 3)
        mesh.AddVertex(-1.0, -1.0, 0.0)
        mesh.AddVertex(1.0, -1.0, 0.0)
        mesh.AddVertex(1.0,  1.0, 0.0)
        mesh.AddVertex(-1.0,  1.0, 0.0)
        mesh.AddVertex(-1.0, -1.0, 1.0)
        mesh.AddVertex(1.0, -1.0, 1.0)
        mesh.AddVertex(1.0,  1.0, 1.0)
        mesh.AddVertex(-1.0,  1.0, 1.0)
        mesh.AddVertex(0.0, -1.0, 0.5)
        mesh.AddVertex(1.0,  0.0, 0.5)
        mesh.AddVertex(0.0,  1.0, 0.5)
        mesh.AddVertex(-1.0,  0.0, 0.5)

        mesh.AddTriangle(0, 1, 8)
        mesh.AddTriangle(1, 5, 8)
        mesh.AddTriangle(5, 4, 8)
        mesh.AddTriangle(4, 0, 8)
        mesh.AddTriangle(1, 2, 9)
        mesh.AddTriangle(2, 6, 9)
        mesh.AddTriangle(6, 5, 9)
        mesh.AddTriangle(5, 1, 9)
        mesh.AddTriangle(2, 3, 10)
        mesh.AddTriangle(3, 7, 10)
        mesh.AddTriangle(7, 6, 10)
        mesh.AddTriangle(6, 2, 10)
        mesh.AddTriangle(3, 0, 11)
        mesh.AddTriangle(0, 4, 11)
        mesh.AddTriangle(4, 7, 11)
        mesh.AddTriangle(7, 3, 11)

        mesh.AddBdrSegment(0, 1, 1)
        mesh.AddBdrSegment(1, 2, 1)
        mesh.AddBdrSegment(2, 3, 1)
        mesh.AddBdrSegment(3, 0, 1)
        mesh.AddBdrSegment(5, 4, 2)
        mesh.AddBdrSegment(6, 5, 2)
        mesh.AddBdrSegment(7, 6, 2)
        mesh.AddBdrSegment(4, 7, 2)
    elif mesh_type == 4:
        mesh = mfem.Mesh(2, 8, 4, 8, 3)

        mesh.AddVertex(-1.0, -1.0, 0.0)
        mesh.AddVertex(1.0, -1.0, 0.0)
        mesh.AddVertex(1.0,  1.0, 0.0)
        mesh.AddVertex(-1.0,  1.0, 0.0)
        mesh.AddVertex(-1.0, -1.0, 1.0)
        mesh.AddVertex(1.0, -1.0, 1.0)
        mesh.AddVertex(1.0,  1.0, 1.0)
        mesh.AddVertex(-1.0,  1.0, 1.0)

        mesh.AddQuad(0, 1, 5, 4)
        mesh.AddQuad(1, 2, 6, 5)
        mesh.AddQuad(2, 3, 7, 6)
        mesh.AddQuad(3, 0, 4, 7)

        mesh.AddBdrSegment(0, 1, 1)
        mesh.AddBdrSegment(1, 2, 1)
        mesh.AddBdrSegment(2, 3, 1)
        mesh.AddBdrSegment(3, 0, 1)
        mesh.AddBdrSegment(5, 4, 2)
        mesh.AddBdrSegment(6, 5, 2)
        mesh.AddBdrSegment(7, 6, 2)
        mesh.AddBdrSegment(4, 7, 2)

    else:
        assert False, "Unrecognized mesh type :" + str(mesh_type)

    mesh.FinalizeTopology()

    return mesh


def run(order=3,
        ref_levels=0,
        static_cond=False,
        mesh_type=4,
        mesh_order=3,
        visualization=True):

    # 2. Construct a quadrilateral or triangular mesh with the topology of a
    #    cylindrical surface.
    mesh = GetMesh(mesh_type)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()

    # 3. Refine the mesh to increase the resolution. In this example we do
    #   'ref_levels' of uniform refinement.
    for l in range(ref_levels):
        mesh.UniformRefinement()

    # 4. Transform the mesh so that it has a more interesting geometry.
    mesh.SetCurvature(mesh_order)

    class cTrans(mfem.VectorPyCoefficient):
        def __init__(self):
            mfem.VectorPyCoefficient.__init__(self, sdim)

        def EvalValue(self, x):
            tol = 1e-6
            theta = 0.0
            if abs(x[1] + 1.0) < tol:
                theta = 0.25 * pi * (x[0] - 2.0)
            elif abs(x[0] - 1.0) < tol:
                theta = 0.25 * pi * x[1]
            elif abs(x[1] - 1.0) < tol:
                theta = 0.25 * pi * (2.0 - x[0])
            elif abs(x[0] + 1.0) < tol:
                theta = 0.25 * pi * (4.0 - x[1])
            else:
                print("side not recognized " +
                      str(x[0]) + " " + str(x[1]) + " " + str(x[2]))

            return (cos(theta), sin(theta),
                    0.25 * (2.0 * x[2] - 1.0) * (cos(theta) + 2.0),)

    trans = cTrans()
    mesh.Transform(trans)

    # 5. Define a finite element space on the mesh. Here we use continuous
    #    Lagrange finite elements of the specified order.
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec)
    print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))

    # 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    #    In this example, the boundary conditions are defined by marking all
    #    the boundary attributes from the mesh as essential (Dirichlet) and
    #    converting them to a list of true dofs.
    ess_tdof_list = mfem.intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr = mfem.intArray(mesh.bdr_attributes.Max())
        ess_bdr.Assign(1)
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 7. Set up the linear form b(.) which corresponds to the right-hand side of
    #    the FEM linear system, which in this case is (1,phi_i) where phi_i are
    #    the basis functions in the finite element fespace.
    b = mfem.LinearForm(fespace)
    one = mfem.ConstantCoefficient(1.0)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
    b.Assemble()

    # 8. Define the solution vector x as a finite element grid function
    #    corresponding to fespace. Initialize x with initial guess of zero,
    #    which satisfies the boundary conditions.
    x = mfem.GridFunction(fespace)
    x.Assign(0.0)

    # 9. Set up the bilinear form a(.,.) on the finite element space
    #    corresponding to the Laplacian operator -Delta, by adding the Diffusion
    #    domain integrator.
    def sigmaFunc(x):
        a = 17.0 - 2.0 * x[0] * (1.0 + x[0])
        s = array([[0.5 + x[0] * x[0] * (8.0 / a - 0.5),
                    x[0] * x[1] * (8.0 / a - 0.5),
                    0.0],
                   [x[0] * x[1] * (8.0 / a - 0.5),
                    0.5 * x[0] * x[0] + 8.0 * x[1] * x[1] / a,
                    0.0],
                   [0.0, 0.0, a/32.]])
        return s

    class cSigma(mfem.MatrixPyCoefficient):
        def __init__(self):
            mfem.MatrixPyCoefficient.__init__(self, sdim)

        def EvalValue(self, x):
            return sigmaFunc(x)

    a = mfem.BilinearForm(fespace)
    sigma = cSigma()
    integ = mfem.DiffusionIntegrator(sigma)
    a.AddDomainIntegrator(integ)

    # 10. Assemble the bilinear form and the corresponding linear system,
    #     applying any necessary transformations such as: eliminating boundary
    #     conditions, applying conforming constraints for non-conforming AMR,
    #     static condensation, etc.
    if static_cond:
        a.EnableStaticCondensation()
    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()

    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

    print("Size of linear system: " + str(A.Height()))

    # 11. Solve the linear system A X = B.
    #     Use a simple symmetric Gauss-Seidel preconditioner with PCG.
    AA = mfem.OperatorHandle2SparseMatrix(A)
    M = mfem.GSSmoother(AA)
    mfem.PCG(A, M, B, X, 1, 200, 1e-12, 0.0)

    # 12. Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, x)

    # 13. Compute error in the solution and its flux
    class cuExact(mfem.PyCoefficient):
        def __init__(self):
            mfem.PyCoefficient.__init__(self)

        def EvalValue(self, x):
            return (0.25 * (2.0 + x[0]) - x[2]) * (x[2] + 0.25 * (2.0 + x[0]))

    def duExactFunc(x):
        return (0.125 * (2.0 + x[0]) * x[1] * x[1],
                -0.125 * (2.0 + x[0]) * x[0] * x[1],
                -2.0 * x[2])

    class cduExact(mfem.VectorPyCoefficient):
        def __init__(self):
            mfem.VectorPyCoefficient.__init__(self, sdim)

        def EvalValue(self, x):
            return duExactFunc(x)

    class cfluxExact(mfem.VectorPyCoefficient):
        def __init__(self):
            mfem.VectorPyCoefficient.__init__(self, sdim)

        def EvalValue(self, x):
            s = sigmaFunc(x)
            du = duExactFunc(x)
            return -s.dot(du)

    uCoef = cuExact()
    err = x.ComputeL2Error(uCoef)

    print("|u - u_h|_2 = " + "{:g}".format(err))

    flux_fespace = mfem.FiniteElementSpace(mesh, fec, 3)
    flux = mfem.GridFunction(flux_fespace)
    x.ComputeFlux(integ, flux)
    flux *= -1.0

    fluxCoef = cfluxExact()
    flux_err = flux.ComputeL2Error(fluxCoef)

    print("|f - f_h|_2 = " + "{:g}".format(flux_err))

    # 14. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g sol.gf".
    mesh.Print("refined.mesh", 8)
    x.Save("sol.gf", 8)

    # 15. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << x
        sol_sock << "window_title 'Solution'\n"
        sol.sock.flush()

        flux_sock = mfem.socketstream("localhost", 19916)
        flux_sock.precision(8)
        flux_sock << "solution\n" << mesh << flux
        flux_sock << "keys vvv\n"
        flux_sock << "window_geometry 402 0 400 350\n"
        flux_sock << "window_title 'Flux'\n"
        flux_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex29 (Solving PDEs on embeded surfaces)')
    parser.add_argument('-mt', '--mesh-type',
                        default=4,
                        action='store', type=int,
                        help="Mesh type: 3 - Triangular, 4 - Quadrilateral.")
    parser.add_argument('-mo', '--mesh-order',
                        default=3,
                        action='store', type=int,
                        help="Geometric order of the curved mesh.")
    parser.add_argument("-r", "--refine",
                        action='store', type=int, default=0,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument('-o', '--order',
                        action='store', default=3, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")
    parser.add_argument('-sc', '--static-condensation',
                        action='store_true',
                        default=False,
                        help="Enable static condensation.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    run(order=args.order,
        ref_levels=args.refine,
        static_cond=args.static_condensation,
        mesh_type=args.mesh_type,
        mesh_order=args.mesh_order,
        visualization=args.visualization,)
