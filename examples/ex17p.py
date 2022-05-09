'''
   MFEM example 17p 

   How to run:
       mpirun -np 4 python <arguments>

   Example of arguments:
       ex17p.py -m beam-tri.mesh
       ex17p.py -m beam-quad.mesh
       ex17p.py -m beam-tet.mesh
       ex17p.py -m beam-hex.mesh
       ex17p.py -m beam-quad.mesh -rs 2 -rp 2 -o 3 -elast
       ex17p.py -m beam-quad.mesh -rs 2 -rp 3 -o 2 -a 1 -k 1
       ex17p.py -m beam-hex.mesh -rs 2 -rp 1 -o 2

'''
import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np


import mfem.par as mfem
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank


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
parser.add_argument('-rs', '--refine-serial',
                    action='store', default=-1, type=int,
                    help="Number of times to refine the mesh uniformly before parallel")
parser.add_argument('-rp', '--refine-parallel',
                    action='store', default=1, type=int,
                    help="Number of times to refine the mesh uniformly after parallel")
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
parser.add_argument('-elast', '--amg-for-elasticity',
                    action='store_true',
                    help='Use the special AMG elasticity solver (GM/LN approaches)',
                    dest='amg_elast', default=False)
parser.add_argument('-sys', '--amg-for-systems',
                    action='store_false',
                    help='Use  standard AMG for systems (unknown approach).',
                    dest='amg_elast', default=True)
parser.add_argument('-vis', '--visualization',
                    action='store_true',
                    help='Enable GLVis visualization')

args = parser.parse_args()
ser_ref_levels = args.refine_serial
par_ref_levels = args.refine_parallel
order = args.order
alpha = args.alpha
kappa = args.kappa
amg_elast = args.amg_elast
visualization = args.visualization
if (kappa < 0):
    kappa = (order+1.)*(order+1.)
    args.kappa = kappa
if (myid == 0):
    parser.print_options(args)
    
device = mfem.Device('cpu')
if myid == 0:            
    device.Print()

# 2. Read the mesh from the given mesh file.
meshfile = expanduser(join(dirname(__file__), '..', 'data', args.mesh))
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
if (mesh.attributes.Max() < 2 or
        mesh.bdr_attributes.Max() < 2):
    if (myid == 0):
        print("\n".join(["Input mesh should have at least two materials and ",
                         "two boundary attributes! (See schematic in ex17.cpp)\n"]))
        sys.exit()

# 3. Refine the mesh to increase the resolution.
if ser_ref_levels < 0:
    ser_ref_levels = int(np.floor(np.log(5000./mesh.GetNE())/np.log(2.)/dim))
for x in range(ser_ref_levels):
    mesh.UniformRefinement()

# Since NURBS meshes do not support DG integrators, we convert them to
# regular polynomial mesh of the specified (solution) order.
if (mesh.NURBSext):
    mesh.SetCurvature(order)
pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
del mesh
for x in range(par_ref_levels):
    pmesh.UniformRefinement()

# 4. Define a DG vector finite element space on the mesh. Here, we use
#    Gauss-Lobatto nodal basis because it gives rise to a sparser matrix
#    compared to the default Gauss-Legendre nodal basis.
fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
fespace = mfem.ParFiniteElementSpace(pmesh, fec, dim, mfem.Ordering.byVDIM)

glob_size = fespace.GlobalTrueVSize()
if (myid == 0):
    print('Number of finite element unknowns: ' + str(glob_size))
    print('Assembling:')

# 5. In this example, the Dirichlet boundary conditions are defined by
#    marking boundary attributes 1 and 2 in the marker Array 'dir_bdr'.
#    These b.c. are imposed weakly, by adding the appropriate boundary
#    integrators over the marked 'dir_bdr' to the bilinear and linear forms.
#    With this DG formulation, there are no essential boundary conditions.
ess_tdof_list = mfem.intArray()
dir_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
dir_bdr.Assign(0)
dir_bdr[0] = 1  # boundary attribute 1 is Dirichlet
dir_bdr[1] = 1  # boundary attribute 2 is Dirichlet

# 6. Define the DG solution vector 'x' as a finite element grid function
#    corresponding to fespace. Initialize 'x' using the 'InitDisplacement'
#    function.
x = mfem.ParGridFunction(fespace)
init_x = InitDisplacement(dim)
x.ProjectCoefficient(init_x)

# 7. Set up the Lame constants for the two materials. They are defined as
#    piece-wise (with respect to the element attributes) constant
#    coefficients, i.e. type PWConstCoefficient.
lamb = mfem.Vector(pmesh.attributes.Max())  # lambda is not possible in python
lamb.Assign(1.0)
lamb[0] = 50.
lambda_c = mfem.PWConstCoefficient(lamb)
mu = mfem.Vector(pmesh.attributes.Max())
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
b = mfem.ParLinearForm(fespace)
if (myid == 0):
    print('r.h.s ...')
integrator = mfem.DGElasticityDirichletLFIntegrator(init_x, lambda_c,
                                                    mu_c, alpha, kappa)
b.AddBdrFaceIntegrator(integrator, dir_bdr)
b.Assemble()

# 9. Set up the bilinear form a(.,.) on the DG finite element space
#    corresponding to the linear elasticity integrator with coefficients
#    lambda and mu as defined above. The additional interior face integrator
#    ensures the weak continuity of the displacement field. The additional
#    boundary face integrator works together with the boundary integrator
#    added to the linear form b(.) to impose weakly the Dirichlet boundary
#    conditions.
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_c, mu_c))

a.AddInteriorFaceIntegrator(
    mfem.DGElasticityIntegrator(lambda_c, mu_c, alpha, kappa))
a.AddBdrFaceIntegrator(
    mfem.DGElasticityIntegrator(lambda_c, mu_c, alpha, kappa), dir_bdr)
if (myid == 0):
    print('matrix ...')
a.Assemble()

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
if (myid == 0):
    print('done.')

# 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
#     solve the system Ax=b with PCG for the symmetric formulation, or GMRES
#     for the non-symmetric.

rtol = 1e-6
amg = mfem.HypreBoomerAMG(A)
if (amg_elast):
    amg.SetElasticityOptions(fespace)
else:
    amg.SetSystemsOptions(dim)

if (alpha == -1.0):
    solver = mfem.CGSolver(A.GetComm())
else:
    solver = mfem.GMRESSolver(A.GetComm())
    solver.SetKDim(50)

solver.SetRelTol(rtol)
solver.SetMaxIter(500)
solver.SetPrintLevel(1)
solver.SetOperator(A)
solver.SetPreconditioner(amg)
solver.Mult(B, X)

# 12. Recover the solution as a finite element grid function 'x'.
a.RecoverFEMSolution(X, b, x)

# 13. Use the DG solution space as the mesh nodal space. This allows us to
#     save the displaced mesh as a curved DG mesh.
pmesh.SetNodalFESpace(fespace)
reference_nodes = mfem.Vector()

if (visualization):
    reference_nodes.Assign(pmesh.GetNodes())
# 14. Save the displaced mesh and minus the solution (which gives the
#     backward displacements to the reference mesh). This output can be
#     viewed later using GLVis
nodes = pmesh.GetNodes()
nodes += x
x.Neg()

smyid = '{:0>6d}'.format(myid)
mesh_name = "mesh."+smyid
sol_name = "sol." + smyid

pmesh.Print(mesh_name, 8)
x.Save(sol_name, 8)

# 15. Visualization: send data by socket to a GLVis server.
if (visualization):
    vis = VisMan("localhost", 19916)
    glvis_keys = "Rjlc" if (dim < 3) else "c"

    vis.NewWindow()
    vis.send_text("parallel " + str(pmesh.GetNRanks()) + " " +
                  str(pmesh.GetMyRank()))
    vis.send_solution(pmesh, x)
    vis.send_text("keys " + glvis_keys)
    vis.send_text("window_title 'Deformed configuration'")
    vis.send_text("plot_caption 'Backward displacement'")
    vis.PositionWindow()
    vis.CloseConnection()

    c = "xyz"
    scalar_dg_space = mfem.ParFiniteElementSpace(pmesh, fec)
    stress = mfem.ParGridFunction(scalar_dg_space)
    stress_c = StressCoefficient(lambda_c, mu_c)

    pmesh.GetNodes().Assign(reference_nodes)
    x.Neg()
    stress_c.SetDisplacement(x)

    def make_plot(si, sj):
        stress_c.SetComponent(si, sj)
        stress.ProjectCoefficient(stress_c)
        vis.NewWindow()
        vis.send_text("parallel " + str(pmesh.GetNRanks()) + " " +
                      str(pmesh.GetMyRank()))
        vis.send_solution(pmesh, stress)
        vis.send_text("keys " + glvis_keys)
        vis.send_text("window_title |Stress" + c[si] + c[sj] + "|")
        vis.PositionWindow()
        vis.CloseConnection()

    for si in range(dim):
        for jj in range(dim-si):
            make_plot(si, si+jj)
