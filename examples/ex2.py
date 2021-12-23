'''
   MFEM example 2

   See c++ version in the MFEM library for more detail
'''
import os
from os.path import expanduser, join, dirname
import numpy as np

import mfem.ser as mfem


def run(order=1, static_cond=False,
        meshfile='', visualization=False,
        device='cpu', pa=False):
    '''
    run ex2
    '''
    device = mfem.Device(device)
    device.Print()

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    #   3. Refine the mesh to increase the resolution. In this example we do
    #      'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
    #      largest number that gives a final mesh with no more than 50,000
    #      elements.

    ref_levels = int(np.floor(np.log(5000. / mesh.GetNE()) / np.log(2.) / dim))
    for x in range(ref_levels):
        mesh.UniformRefinement()

    # 5. Define a finite element space on the mesh. Here we use vector finite
    #   elements, i.e. dim copies of a scalar finite element space. The vector
    #   dimension is specified by the last argument of the FiniteElementSpace
    #   constructor. For NURBS meshes, we use the (degree elevated) NURBS space
    #   associated with the mesh nodes.
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec, dim)

    print("Number of finite element unknowns: " + str(fespace.GetTrueVSize()))

    # 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    #    In this example, the boundary conditions are defined by marking only
    #    boundary attribute 1 from the mesh as essential and converting it to a
    #    list of true dofs.
    ess_tdof_list = mfem.intArray()
    ess_bdr = mfem.intArray([1] + [0] * (mesh.bdr_attributes.Max() - 1))
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)
    # 7. Set up the linear form b(.) which corresponds to the right-hand side of
    #    the FEM linear system. In this case, b_i equals the boundary integral
    #    of f*phi_i where f represents a "pull down" force on the Neumann part
    #    of the boundary and phi_i are the basis functions in the finite element
    #    fespace. The force is defined by the VectorArrayCoefficient object f,
    #    which is a vector of Coefficient objects. The fact that f is non-zero
    #    on boundary attribute 2 is indicated by the use of piece-wise constants
    #    coefficient for its last component.
    f = mfem.VectorArrayCoefficient(dim)
    for i in range(dim - 1):
        f.Set(i, mfem.ConstantCoefficient(0.0))

    pull_force = mfem.Vector([0] * mesh.bdr_attributes.Max())
    pull_force[1] = -1.0e-2
    f.Set(dim - 1, mfem.PWConstCoefficient(pull_force))

    b = mfem.LinearForm(fespace)
    b.AddBoundaryIntegrator(mfem.VectorBoundaryLFIntegrator(f))
    print('r.h.s...')
    b.Assemble()
    # 8. Define the solution vector x as a finite element grid function
    #    corresponding to fespace. Initialize x with initial guess of zero,
    #    which satisfies the boundary conditions.
    x = mfem.GridFunction(fespace)
    x.Assign(0.0)
    # 9. Set up the bilinear form a(.,.) on the finite element space
    #    corresponding to the linear elasticity integrator with piece-wise
    #    constants coefficient lambda and mu.
    #    (Here lambda is renamed by lamb since it is reserved in python)
    lamb = mfem.Vector(mesh.attributes.Max())
    lamb.Assign(1.0)
    lamb[0] = lamb[1] * 50
    lambda_func = mfem.PWConstCoefficient(lamb)
    mu = mfem.Vector(mesh.attributes.Max())
    mu.Assign(1.0)
    mu[0] = mu[1] * 50
    mu_func = mfem.PWConstCoefficient(mu)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_func, mu_func))
    # 10. Assemble the bilinear form and the corresponding linear system,
    #     applying any necessary transformations such as: eliminating boundary
    #     conditions, applying conforming constraints for non-conforming AMR,
    #     static condensation, etc.
    print('matrix...')
    if (static_cond):
        a.EnableStaticCondensation()
    a.Assemble()

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
    print('...done')  # Here, original version calls hegith, which is not
    # defined in the header...!?
    print("Size of linear system: " + str(A.Height()))

    # 10. Solve
    AA = mfem.OperatorHandle2SparseMatrix(A)
    M = mfem.GSSmoother(AA)
    mfem.PCG(AA, M, B, X, 1, 500, 1e-8, 0.0)

    # 11. Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, x)

    # 13. For non-NURBS meshes, make the mesh curved based on the finite element
    #     space. This means that we define the mesh elements through a fespace
    #     based transformation of the reference element. This allows us to save
    #     the displaced mesh as a curved mesh when using high-order finite
    #     element displacement field. We assume that the initial mesh (read from
    #     the file) is not higher order curved mesh compared to the chosen FE
    #     space.
    if not mesh.NURBSext:
        mesh.SetNodalFESpace(fespace)

    # 14. Save the displaced mesh and the inverted solution (which gives the
    #     backward displacements to the original grid). This output can be
    #     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
    nodes = mesh.GetNodes()
    nodes += x
    x *= -1
    mesh.Print('displaced.mesh', 8)
    x.Save('sol.gf', 8)

    if (visualization):
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock.send_solution(mesh, x)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser
    parser = ArgParser(description='Ex2 (Linear elasticity)')
    parser.add_argument('-m', '--mesh',
                        default='beam-tri.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
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

    args = parser.parse_args()
    parser.print_options(args)

    order = args.order
    static_cond = args.static_condensation

    meshfile = expanduser(
        join(dirname(__file__), '..', 'data', args.mesh))

    visualization = args.visualization
    device = args.device
    pa = args.partial_assembly

    run(order=order,
        static_cond=static_cond,
        meshfile=meshfile,
        visualization=visualization,
        device=device,
        pa=pa)
