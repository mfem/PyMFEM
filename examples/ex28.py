'''
   MFEM example 28
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.ser as mfem
from mfem.ser import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array, floor, log


def build_trapezoid_mesh(offset):
    assert offset < 0.9, "offset is too large!"

    dimension = 2
    nvt = 4  # vertices
    nbe = 4  # num boundary elements
    mesh = mfem.Mesh(dimension, nvt, 1, nbe)

    mesh.AddVertex(0.0, 0.0)
    mesh.AddVertex(1.0, 0.0)
    mesh.AddVertex(offset, 1.0)
    mesh.AddVertex(1.0, 1.0)

    #  element
    mesh.AddQuad([0, 1, 3, 2], 1)

    #  boundary
    mesh.AddBdrSegment([0, 1], 1)
    mesh.AddBdrSegment([1, 3], 2)
    mesh.AddBdrSegment([2, 3], 3)
    mesh.AddBdrSegment([0, 2], 4)

    mesh.FinalizeQuadMesh(1, 0, True)

    return mesh


def run(order=1,
        offset=0.3,
        visit=False,
        visualization=True):

    # 2. Build a trapezoidal mesh with a single quadrilateral element, where
    #    'offset' determines how far off it is from a rectangle.
    mesh = build_trapezoid_mesh(offset)
    dim = mesh.Dimension()

    # 3. Refine the mesh to increase the resolution. In this example we do
    #    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    #    largest number that gives a final mesh with no more than 1,000
    #    elements.
    ref_levels = int(floor(log(1000./mesh.GetNE())/log(2.)/dim))
    for l in range(ref_levels):
        mesh.UniformRefinement()

    # 4. Define a finite element space on the mesh. Here we use vector finite
    #    elements, i.e. dim copies of a scalar finite element space. The vector
    #    dimension is specified by the last argument of the FiniteElementSpace
    #    constructor.

    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec, dim)
    print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))
    print("Assembling matrix and r.h.s... ")

    # 5. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, there are no essential boundary
    #    conditions in the usual sense, but we leave the machinery here for
    #    users to modify if they wish.
    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray(mesh.bdr_attributes.Max())
    ess_bdr.Assign(0)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 6. Set up the linear form b(.) which corresponds to the right-hand side of
    #    the FEM linear system. In this case, b_i equals the boundary integral
    #    of f*phi_i where f represents a "push" force on the right side of the
    #    trapezoid.
    f = mfem.VectorArrayCoefficient(dim)
    for i in range(dim-1):
        f.Set(i, mfem.ConstantCoefficient(0.0))
    push_force = mfem.Vector(mesh.bdr_attributes.Max())
    push_force.Assign(0.0)
    push_force[1] = -5.0e-2
    f.Set(0, mfem.PWConstCoefficient(push_force))

    b = mfem.LinearForm(fespace)
    b.AddBoundaryIntegrator(mfem.VectorBoundaryLFIntegrator(f))
    b.Assemble()

    # 7. Define the solution vector x as a finite element grid function
    #    corresponding to fespace.
    x = mfem.GridFunction(fespace)
    x.Assign(0.0)

    # 8. Set up the bilinear form a(.,.) on the finite element space
    #    corresponding to the linear elasticity integrator with piece-wise
    #    constants coefficient lambda and mu. We use constant coefficients,
    #    but see ex2 for how to set up piecewise constant coefficients based
    #    on attribute.
    llambda = mfem.Vector(mesh.attributes.Max())
    llambda.Assign(1.0)
    lambda_func = mfem.PWConstCoefficient(llambda)
    mu = mfem.Vector(mesh.attributes.Max())
    mu.Assign(1.0)
    mu_func = mfem.PWConstCoefficient(mu)

    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_func, mu_func))

    # 9. Assemble the bilinear form and the corresponding linear system,
    #    applying any necessary transformations such as: eliminating boundary
    #    conditions, applying conforming constraints for non-conforming AMR,
    #    static condensation, etc.
    a.Assemble()

    A = mfem.SparseMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
    print("done.")
    print("Size of linear system: " + str(A.Height()))

    # 10. Set up constraint matrix to constrain normal displacement (but
    #     allow tangential displacement) on specified boundaries.
    constraint_atts = mfem.intArray([1, 4])
    lagrange_rowstarts = mfem.intArray()
    local_constraints = mfem.BuildNormalConstraints(fespace,
                                                    constraint_atts,
                                                    lagrange_rowstarts)

    # 11. Define and apply an iterative solver for the constrained system
    #     in saddle-point form with a Gauss-Seidel smoother for the
    #     displacement block.
    M = mfem.GSSmoother(A)
    solver = mfem.SchurConstrainedSolver(A, local_constraints, M)
    solver.SetRelTol(1e-5)
    solver.SetMaxIter(2000)
    solver.SetPrintLevel(1)
    solver.Mult(B, X)

    # 12. Recover the solution as a finite element grid function. Move the
    #     mesh to reflect the displacement of the elastic body being
    #     simulated, for purposes of output.
    a.RecoverFEMSolution(X, b, x)
    mesh.SetNodalFESpace(fespace)
    nodes = mesh.GetNodes()
    nodes += x

    #  13. Save the refined mesh and the solution in VisIt format.
    if visit:
        visit_dc = mfem.VisItDataCollection("ex28", mesh)
        visit_dc.SetLevelsOfDetail(4)
        visit_dc.RegisterField("displacement", x)
        visit_dc.Save()

    # 14. Save the displaced mesh and the inverted solution (which gives the
    #     backward displacements to the original grid). This output can be
    #     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
    x *= -1  # sign convention for GLVis displacements
    mesh.Print("displaced.mesh", 8)
    x.Save("sol.gf", 8)

    # 15. Send the above data by socket to a GLVis server.  Use the "n" and "b"
    #     keys in GLVis to visualize the displacements.
    if visualization:
        sol_sock = mfem.socketstream('localhost', 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << x
        sol_sock.flush


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex28 (Constraints and sliding boundary conditions')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-offset', '--offset',
                        action='store', default=0.3, type=float,
                        help="How much to offset the trapezoid.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-visit', '--visit-datafiles',
                        action='store_true',
                        help="Save data files for VisIt (visit.llnl.gov) visualization.")

    args = parser.parse_args()
    parser.print_options(args)

    run(order=args.order,
        offset=args.offset,
        visit=args.visit_datafiles,
        visualization=args.visualization,)
