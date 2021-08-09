'''
   MFEM example 17 

   How to run:
       python <arguments>

   Example of arguments:
       ex17.py -m beam-tri.mesh
       ex17.py -m beam-quad.mesh
       ex17.py -m beam-tet.mesh
       ex17.py -m beam-hex.mesh
       ex17.py -m beam-quad.mesh -r 2 -o 3
       ex17.py -m beam-quad.mesh -r 2 -o 2 -a 1 -k 1
       ex17.py -m beam-hex.mesh -r 2 -o 2

'''
import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from mfem import path

import mfem.ser as mfem
from mfem.ser import intArray

#


class InitDisplacement(mfem.VectorPyCoefficient):
    def __init__(self, dim):
        self.dim = dim
        mfem.VectorPyCoefficient.__init__(self, dim)

    def EvalValue(self, x):
        u = [0.0]*dim
        u[-1] = -0.2*x[0]
        return tuple(u)


class StressCoefficient(mfem.PyCoefficientBase):
    def __init__(self, lambda_, mu_, si=0, sj=0):
        super(StressCoefficient, self).__init__(0)
        self.lam = lambda_   # coefficient
        self.mu = mu_       # coefficient
        self.si = si
        self.sj = sj     # component
        self.u = None   # displacement GridFunction
        self.grad = mfem.DenseMatrix()

    def SetComponent(self, i, j):
        self.si = i
        self.sj = j

    def SetDisplacement(self, u):
        self.u = u

    def Eval(self, T, ip):
        si, sj = self.si, self.sj
        L = self.lam.Eval(T, ip)
        M = self.mu.Eval(T, ip)
        self.u.GetVectorGradient(T, self.grad)
        if (self.si == self.sj):
            div_u = self.grad.Trace()
            return L * div_u + 2 * M * self.grad[si, si]
        else:
            return M * (self.grad[si, sj] + self.grad[sj, si])


class VisMan(object):
    def __init__(self, vishost, visport):
        self.host = vishost
        self.port = visport
        self.socks = []
        self.output = None
        self.win_x = 0
        self.win_y = 0
        self.win_w = 200  # window width
        self.win_h = 150  # window height
        self.stride_x = self.win_w
        self.stride_y = self.win_h + 20
        self.win_nx = 4  # number of windows in a row
        self.sid = 0

    def NewWindow(self):
        self.socks.append(mfem.socketstream(self.host, self.port))
        self.output = self.socks[-1]
        self.output.precision(8)
        self.socks
        self.sid = self.sid + 1

    def CloseConnection(self):
        self.socks = []
        del self.output
        self.output = None

    def PositionWindow(self):
        if self.output is None:
            return

        sid = self.sid
        command = ("window_geometry " +
                   str(self.win_x + self.stride_x*(sid % self.win_nx)) +
                   ' ' +
                   str(self.win_y + self.stride_y*(sid/self.win_nx)) +
                   ' ' + str(self.win_w) + ' ' + str(self.win_h))
        self.output.send_text(command)
        self.output.flush()

    def send_solution(self, mesh, x):
        if self.output is None:
            return
        self.output.send_solution(mesh, x)

    def send_text(self, x):
        if self.output is None:
            return
        self.output.send_text(x)

    def flush(self):
        if self.output is None:
            return
        self.output.flush()


parser = ArgParser(description='Ex17')
parser.add_argument('-m', '--mesh',
                    default='beam-tri.mesh',
                    action='store', type=str,
                    help='Mesh file to use.')
parser.add_argument('-r', '--refine',
                    action='store', default=-1, type=int,
                    help="Number of times to refine the mesh uniformly, -1 for auto.")
parser.add_argument('-o', '--order',
                    action='store', default=1, type=int,
                    help="Finite element order (polynomial degree)")
parser.add_argument('-a', '--alpha',
                    action='store', default=-1.0, type=float,
                    help='\n'.join(["One of the two DG penalty parameters, typically +1/-1."
                                    " See the documentation of class DGElasticityIntegrator."]))
parser.add_argument('-k', '--kappa',
                    action='store', default=-1.0, type=float,
                    help='\n'.join(["One of the two DG penalty parameters, should be positve."
                                    " Negative values are replaced with (order+1)^2."]))
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')

args = parser.parse_args()
ref_levels = args.refine
order = args.order
alpha = args.alpha
kappa = args.kappa
visualization = args.visualization
if (kappa < 0):
    kappa = (order+1.)*(order+1.)
    args.kappa = kappa
parser.print_options(args)
# 2. Read the mesh from the given mesh file
meshfile = expanduser(
    join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
if (mesh.attributes.Max() < 2 or
        mesh.bdr_attributes.Max() < 2):
    print("\n".join(["Input mesh should have at least two materials and ",
                     "two boundary attributes! (See schematic in ex17.cpp)\n"]))
    sys.exit()

# 3. Refine the mesh to increase the resolution.
ref_levels = int(np.floor(np.log(5000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ref_levels):
    mesh.UniformRefinement()

# Since NURBS meshes do not support DG integrators, we convert them to
# regular polynomial mesh of the specified (solution) order.
if (mesh.NURBSext):
    mesh.SetCurvature(order)

# 4. Define a DG vector finite element space on the mesh. Here, we use
#    Gauss-Lobatto nodal basis because it gives rise to a sparser matrix
#    compared to the default Gauss-Legendre nodal basis.
fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
fespace = mfem.FiniteElementSpace(mesh, fec, dim)
print('Number of finite element unknowns: ' + str(fespace.GetVSize()))
print('Assembling:')

# 5. In this example, the Dirichlet boundary conditions are defined by
#    marking boundary attributes 1 and 2 in the marker Array 'dir_bdr'.
#    These b.c. are imposed weakly, by adding the appropriate boundary
#    integrators over the marked 'dir_bdr' to the bilinear and linear forms.
#    With this DG formulation, there are no essential boundary conditions.
ess_tdof_list = intArray()
dir_bdr = intArray(mesh.bdr_attributes.Max())
dir_bdr.Assign(0)
dir_bdr[0] = 1  # boundary attribute 1 is Dirichlet
dir_bdr[1] = 1  # boundary attribute 2 is Dirichlet

# 6. Define the DG solution vector 'x' as a finite element grid function
#    corresponding to fespace. Initialize 'x' using the 'InitDisplacement'
#    function.
x = mfem.GridFunction(fespace)
init_x = InitDisplacement(dim)
x.ProjectCoefficient(init_x)

# 7. Set up the Lame constants for the two materials. They are defined as
#    piece-wise (with respect to the element attributes) constant
#    coefficients, i.e. type PWConstCoefficient.
lamb = mfem.Vector(mesh.attributes.Max())  # lambda is not possible in python
lamb.Assign(1.0)
lamb[0] = 50.
lambda_c = mfem.PWConstCoefficient(lamb)
mu = mfem.Vector(mesh.attributes.Max())
mu.Assign(1.0)
mu[0] = 50.0
mu_c = mfem.PWConstCoefficient(mu)

# 8. Set up the linear form b(.) which corresponds to the right-hand side of
#    the FEM linear system. In this example, the linear form b(.) consists
#    only of the terms responsible for imposing weakly the Dirichlet
#    boundary conditions, over the attributes marked in 'dir_bdr'. The
#    values for the Dirichlet boundary condition are taken from the
#    VectorFunctionCoefficient 'x_init' which in turn is based on the
#    function 'InitDisplacement'.
b = mfem.LinearForm(fespace)
print('r.h.s ...')
integrator = mfem.DGElasticityDirichletLFIntegrator(
    init_x, lambda_c, mu_c, alpha, kappa)
b.AddBdrFaceIntegrator(integrator, dir_bdr)
b.Assemble()

# 9. Set up the bilinear form a(.,.) on the DG finite element space
#    corresponding to the linear elasticity integrator with coefficients
#    lambda and mu as defined above. The additional interior face integrator
#    ensures the weak continuity of the displacement field. The additional
#    boundary face integrator works together with the boundary integrator
#    added to the linear form b(.) to impose weakly the Dirichlet boundary
#    conditions.
a = mfem.BilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_c, mu_c))

a.AddInteriorFaceIntegrator(
    mfem.DGElasticityIntegrator(lambda_c, mu_c, alpha, kappa))
a.AddBdrFaceIntegrator(mfem.DGElasticityIntegrator(
    lambda_c, mu_c, alpha, kappa), dir_bdr)
print('matrix ...')
a.Assemble()

A = mfem.SparseMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
print('...done')


A.PrintInfo()
'''
   Note: extension of ostream &

   A.PrintInfo()               # output to std::cout
   A.PrintInfo("matrix info")  # output to file 

   Above two are the same as

   from mfem._ser.io_stream import STDOUT, wFILE
   A.PrintInfo(STDOUT)
   A.PrintInfo(wFILE("matrix_info"))

'''

# 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
#     solve the system Ax=b with PCG for the symmetric formulation, or GMRES
#     for the non-symmetric.
M = mfem.GSSmoother(A)
rtol = 1e-6
if (alpha == -1.0):
    mfem.PCG(A, M, B, X, 3, 5000, rtol*rtol, 0.0)
else:
    mfem.GMRES(A, M, B, X, 3, 5000, 50, rtol*rtol, 0.0)

# 12. Recover the solution as a finite element grid function 'x'.
a.RecoverFEMSolution(X, b, x)

# 13. Use the DG solution space as the mesh nodal space. This allows us to
#     save the displaced mesh as a curved DG mesh.
mesh.SetNodalFESpace(fespace)
reference_nodes = mfem.Vector()
if (visualization):
    reference_nodes.Assign(mesh.GetNodes())
# 14. Save the displaced mesh and minus the solution (which gives the
#     backward displacements to the reference mesh). This output can be
#     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
nodes = mesh.GetNodes()
nodes += x
x.Neg()

mesh.Print('displaced.mesh', 8)
x.Save('sol.gf', 8)

# 15. Visualization: send data by socket to a GLVis server.
if (visualization):
    vis = VisMan("localhost", 19916)
    glvis_keys = "Rjlc" if (dim < 3) else "c"

    vis.NewWindow()
    vis.send_solution(mesh, x)
    vis.send_text("keys " + glvis_keys)
    vis.send_text("window_title 'Deformed configuration'")
    vis.send_text("plot_caption 'Backward displacement'")
    vis.PositionWindow()
    vis.CloseConnection()

    c = "xyz"
    scalar_dg_space = mfem.FiniteElementSpace(mesh, fec)
    stress = mfem.GridFunction(scalar_dg_space)
    stress_c = StressCoefficient(lambda_c, mu_c)

    mesh.GetNodes().Assign(reference_nodes)
    x.Neg()
    stress_c.SetDisplacement(x)

    def make_plot(si, sj):
        stress_c.SetComponent(si, sj)
        stress.ProjectCoefficient(stress_c)
        vis.NewWindow()
        vis.send_solution(mesh, stress)
        vis.send_text("keys " + glvis_keys)
        vis.send_text("window_title |Stress" + c[si] + c[sj] + "|")
        vis.PositionWindow()
        vis.CloseConnection()

    for si in range(dim):
        for jj in range(dim-si):
            make_plot(si, si+jj)
