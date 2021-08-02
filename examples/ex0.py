'''
   MFEM example 0 (converted from ex0.cpp)

   See c++ version in the MFEM library for more detail

   How to run:
      python <arguments>

   Example of arguments:
      ex1.py -m star.mesh
      ex1.py -m fichera.mesh -o 2

   Description: This example code demonstrates the most basic usage of MFEM to
                define a simple finite element discretization of the Laplace
                problem -Delta u = 1 with zero Dirichlet boundary conditions.
                General 2D/3D mesh files and finite element polynomial degrees
                can be specified by command line options.

'''
import os
from os.path import expanduser, join
import numpy as np

import mfem.ser as mfem


def run(order=1,  meshfile=''):
    '''
    run ex0
    '''

    #  2. Read the mesh from the given mesh file, and refine once uniformly.
    mesh = mfem.Mesh(meshfile, 1, 1)
    mesh.UniformRefinement()

    # 3. Define a finite element space on the mesh. Here we use H1 continuous
    #    high-order Lagrange finite elements of the given order.
    fec = mfem.H1_FECollection(order,  mesh.Dimension())
    fespace = mfem.FiniteElementSpace(mesh, fec)
    print('Number of finite element unknowns: ' +
          str(fespace.GetTrueVSize()))

    # 4. Extract the list of all the boundary DOFs. These will be marked as
    #    Dirichlet in order to enforce zero boundary conditions.
    boundary_dofs = mfem.intArray()
    fespace.GetBoundaryTrueDofs(boundary_dofs)

    # 5. Define the solution x as a finite element grid function in fespace. Set
    #    the initial guess to zero, which also sets the boundary conditions.
    x = mfem.GridFunction(fespace)
    x.Assign(0.0)

    # 6. Set up the linear form b(.) corresponding to the right-hand side.
    one = mfem.ConstantCoefficient(1.0)
    b = mfem.LinearForm(fespace)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
    b.Assemble()

    # 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
    a.Assemble()

    # 8. Form the linear system A X = B. This includes eliminating boundary
    #    conditions, applying AMR constraints, and other transformations.
    A = mfem.SparseMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(boundary_dofs, x, b, A, X, B)
    print("Size of linear system: " + str(A.Height()))

    # 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
    M = mfem.GSSmoother(A)
    mfem.PCG(A, M, B, X, 1, 200, 1e-12, 0.0)

    # 10. Recover the solution x as a grid function and save to file. The output
    #     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
    a.RecoverFEMSolution(X, b, x)
    x.Save('sol.gf')
    mesh.Save('mesh.mesh')


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex1 (Laplace Problem)')
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")

    args = parser.parse_args()
    parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    run(order=order,
        meshfile=meshfile)
