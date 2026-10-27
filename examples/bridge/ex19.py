"""Standalone serial MFEM example 19 with native NSB callbacks."""

from pathlib import Path
from types import SimpleNamespace
from numba import njit, types
import nsb_rt as nsb
import mfem.ser as mfem
from mfem.ser import intArray
from mfem.common.arg_parser import ArgParser

registration = mfem.get_bridge_registration()


@njit
def block_view(vector, start, end):
    data = vector.GetDataArray()
    assert 0 <= start <= end <= data.size
    return mfem.Vector(nsb.as_cpointer(data[start:end]), end - start)


@njit
def apply_preconditioner(k, y, offsets, mass_pcg, stiff_pcg, jacobian,
                         gamma, scratch, scratch2):
    # Return normally from the allocating frame even if a native call fails,
    # so Numba releases the small C++ view wrappers before failure propagates.
    try:
        offset_data = offsets.GetDataArray()
        assert offset_data.size >= 3
        first, split, end = offset_data[0], offset_data[1], offset_data[2]
        # Check every bound before allocating frame-owned view wrappers.
        assert 0 <= first <= split <= end <= k.Size()
        assert end <= y.Size()
        disp_in = block_view(k, first, split)
        pres_in = block_view(k, split, end)
        disp_out = block_view(y, first, split)
        pres_out = block_view(y, split, end)
        mass_pcg.Mult(pres_in, pres_out)
        pressure = pres_out.GetDataArray()
        for i in range(pressure.size):
            pressure[i] *= -gamma
        jacobian.GetBlock(0, 1).Mult(pres_out, scratch)
        displacement = disp_in.GetDataArray()
        temp = scratch.GetDataArray()
        temp2 = scratch2.GetDataArray()
        for i in range(temp2.size):
            temp2[i] = displacement[i] - temp[i]
        stiff_pcg.Mult(scratch2, disp_out)
        return True
    except Exception:
        return False


def get_directors():
    @nsb.director(state={"point": mfem.Vector}, fallback="silent")
    class InitialDeformation(mfem.VectorCoefficient):
        def __init__(self, dim):
            super().__init__(dim)
            self.point = mfem.Vector(dim)

        @nsb.override
        def Eval(self, y, T, ip):
            T.Transform(ip, self.point)
            x = self.point.GetDataArray()
            y.SetSize(x.size)
            out = y.GetDataArray()
            for i in range(x.size):
                out[i] = x[i]
            out[1] += 0.25 * x[0]

    @nsb.director(state={"point": mfem.Vector}, fallback="silent")
    class ReferenceConfiguration(mfem.VectorCoefficient):
        def __init__(self, dim):
            super().__init__(dim)
            self.point = mfem.Vector(dim)

        @nsb.override
        def Eval(self, y, T, ip):
            T.Transform(ip, self.point)
            x = self.point.GetDataArray()
            y.SetSize(x.size)
            out = y.GetDataArray()
            for i in range(x.size):
                out[i] = x[i]

    @nsb.director(state={"jacobian": mfem.BlockOperator,
                         "block_offsets": mfem.intArray, "gamma": types.float64,
                         "mass_pcg": mfem.CGSolver, "stiff_pcg": mfem.GMRESSolver,
                         "temp": mfem.Vector, "temp2": mfem.Vector}, fallback="silent")
    class JacobianPreconditioner(mfem.Solver):
        def __init__(self, spaces, mass, offsets):
            super().__init__(offsets[2])
            self.pressure_mass = mass
            self.block_offsets = offsets
            self.gamma = 0.00001
            self.temp = mfem.Vector(offsets[1] - offsets[0])
            self.temp2 = mfem.Vector(offsets[1] - offsets[0])
            self.mass_prec = mfem.GSSmoother(mass)
            self.mass_pcg = mfem.CGSolver()
            self.mass_pcg.SetRelTol(1e-12)
            self.mass_pcg.SetAbsTol(1e-12)
            self.mass_pcg.SetMaxIter(200)
            self.mass_pcg.SetPrintLevel(0)
            self.mass_pcg.SetPreconditioner(self.mass_prec)
            self.mass_pcg.SetOperator(mass)
            self.mass_pcg.iterative_mode = False
            self.stiff_prec = mfem.GSSmoother()
            self.stiff_pcg = mfem.GMRESSolver()
            self.stiff_pcg.SetRelTol(1e-8)
            self.stiff_pcg.SetAbsTol(1e-8)
            self.stiff_pcg.SetMaxIter(200)
            self.stiff_pcg.SetPrintLevel(0)
            self.stiff_pcg.SetPreconditioner(self.stiff_prec)
            self.stiff_pcg.iterative_mode = False

        @nsb.override
        def SetOperator(self, op):
            # Hform owns this Jacobian. Every rebuild must be followed by this
            # update before Mult uses its blocks or the stiffness solver.
            self.jacobian = nsb.dynamic_cast(op, mfem.BlockOperator)
            self.stiff_pcg.SetOperator(self.jacobian.GetBlock(0, 0))

        @nsb.override
        def Mult(self, k, y):
            success = apply_preconditioner(k, y, self.block_offsets,
                                          self.mass_pcg, self.stiff_pcg,
                                          self.jacobian, self.gamma,
                                          self.temp, self.temp2)
            assert success, "Jacobian preconditioner failed"

    @nsb.director(state={"Hform": mfem.BlockNonlinearForm}, fallback="silent")
    class RubberOperator(mfem.Operator):
        def __init__(self, spaces, ess_bdr, block_offsets,
                     rel_tol, abs_tol, iter, mu):
            super().__init__(spaces[0].GetVSize() + spaces[1].GetVSize())
            self.spaces = spaces
            self.mu = mfem.ConstantCoefficient(mu)
            self.block_offsets = block_offsets
            self.Hform = mfem.BlockNonlinearForm(spaces)
            self.Hform.AddDomainIntegrator(mfem.IncompressibleNeoHookeanIntegrator(self.mu))
            self.Hform.SetEssentialBC(ess_bdr, (None, None))
            a = mfem.BilinearForm(spaces[1])
            one = mfem.ConstantCoefficient(1.0)
            a.AddDomainIntegrator(mfem.MassIntegrator(one))
            a.Assemble()
            a.Finalize()
            # Retain the form as the single matrix owner. SpMat() is borrowed;
            # Python setup need not transfer ownership through LoseMat().
            self.pressure_mass_form = a
            self.pressure_mass_coefficient = one
            pressure_mass = a.SpMat()
            self.j_prec = JacobianPreconditioner(spaces, pressure_mass, block_offsets)
            self.j_solver = mfem.GMRESSolver()
            self.j_solver.iterative_mode = False
            self.j_solver.SetRelTol(1e-12)
            self.j_solver.SetAbsTol(1e-12)
            self.j_solver.SetMaxIter(300)
            self.j_solver.SetPrintLevel(-1)
            self.j_monitor = GeneralResidualMonitor("  GMRES", 3)
            self.j_solver.SetMonitor(self.j_monitor)
            self.j_solver.SetPreconditioner(self.j_prec)
            self.newton_solver = mfem.NewtonSolver()
            self.newton_solver.iterative_mode = True
            self.newton_solver.SetSolver(self.j_solver)
            self.newton_solver.SetOperator(self)
            self.newton_solver.SetPrintLevel(-1)
            self.newton_solver.SetRelTol(rel_tol)
            self.newton_solver.SetAbsTol(abs_tol)
            self.newton_solver.SetMaxIter(iter)
            self.n_monitor = GeneralResidualMonitor("Newton", 1)
            self.newton_solver.SetMonitor(self.n_monitor)

        def Solve(self, xp):
            self.newton_solver.Mult(mfem.Vector(), xp)
            assert self.newton_solver.GetConverged(), "Newton Solver did not converge."

        @nsb.override
        def Mult(self, k, y):
            self.Hform.Mult(k, y)

        @nsb.override
        def GetGradient(self, xp):
            # Borrowed native return; Hform stays alive through the Newton solve.
            # A subsequent rebuild invalidates references to the old gradient.
            return self.Hform.GetGradient(xp)

        def shutdown(self):
            # Release native consumers before detaching callbacks and their state.
            self.newton_solver = None
            self.j_solver = None
            self.j_prec.close()
            self.close()

    return InitialDeformation, ReferenceConfiguration, JacobianPreconditioner, RubberOperator


class GeneralResidualMonitor(mfem.IterativeSolverController):
    def __init__(self, prefix, print_level):
        mfem.IterativeSolverController.__init__(self)
        self.print_level = print_level
        self.prefix = prefix
        self.norm = 0

    def MonitorResidual(self, it, norm, r, final):
        txt = ''
        if (self.print_level == 1 or
                (self.print_level == 3 and (final or it == 0))):
            txt = (txt + self.prefix + " iteration " + '{:>2}'.format(it) +
                   " : ||r|| = " + "{:g}".format(norm))
            if it > 0:
                txt = txt + ",  ||r||/||r_0|| = " + \
                    "{:g}".format(norm/self.norm0)
            else:
                self.norm0 = norm
            print(txt)

def ex19_main(args):
    ref_levels = args.refine
    order = args.order
    visualization = args.visualization
    mu = args.shear_modulus
    newton_rel_tol = args.relative_tolerance
    newton_abs_tol = args.absolute_tolerance
    newton_iter = args.newton_iterations

    if order < 2 or ref_levels < 0 or mu <= 0:
        raise ValueError("ex19 requires order >= 2, nonnegative refinement, and mu > 0")


    meshfile = Path(args.mesh).expanduser()
    if not meshfile.is_file():
        meshfile = Path(__file__).resolve().parents[2] / 'data' / meshfile
    meshfile = str(meshfile)
    mesh = mfem.Mesh(meshfile, 1, 1)
    InitialDeformation, ReferenceConfiguration, _, RubberOperator = get_directors()
    dim = mesh.Dimension()
    if dim not in (2, 3) or mesh.SpaceDimension() != dim:
        raise ValueError("ex19 requires a full-dimensional 2D or 3D mesh")

    for lev in range(ref_levels):
        mesh.UniformRefinement()

    #  4. Define the shear modulus for the incompressible Neo-Hookean material
    c_mu = mfem.ConstantCoefficient(mu)

    #  5. Define the finite element spaces for displacement and pressure
    #     (Taylor-Hood elements). By default, the displacement (u/x) is a second
    #     order vector field, while the pressure (p) is a linear scalar function.
    quad_coll = mfem.H1_FECollection(order, dim)
    lin_coll = mfem.H1_FECollection(order-1, dim)

    R_space = mfem.FiniteElementSpace(
        mesh, quad_coll, dim, mfem.Ordering.byVDIM)
    W_space = mfem.FiniteElementSpace(mesh, lin_coll)

    spaces = [R_space, W_space]
    R_size = R_space.GetVSize()
    W_size = W_space.GetVSize()

    #   6. Define the Dirichlet conditions (set to boundary attribute 1 and 2)
    ess_bdr_u = mfem.intArray(R_space.GetMesh().bdr_attributes.Max())
    ess_bdr_p = mfem.intArray(W_space.GetMesh().bdr_attributes.Max())
    ess_bdr_u.Assign(0)
    ess_bdr_u[0] = 1
    ess_bdr_u[1] = 1
    ess_bdr_p.Assign(0)
    ess_bdr = [ess_bdr_u, ess_bdr_p]

    print("***********************************************************")
    print("dim(u) = " + str(R_size))
    print("dim(p) = " + str(W_size))
    print("dim(u+p) = " + str(R_size + W_size))
    print("***********************************************************")

    block_offsets = intArray([0, R_size, W_size])
    block_offsets.PartialSum()

    xp = mfem.BlockVector(block_offsets)

    #  9. Define grid functions for the current configuration, reference
    #     configuration, final deformation, and pressure
    x_gf = mfem.GridFunction(R_space)
    x_ref = mfem.GridFunction(R_space)
    x_def = mfem.GridFunction(R_space)
    p_gf = mfem.GridFunction(W_space)

    x_gf.MakeRef(R_space, xp.GetBlock(0), 0)
    p_gf.MakeRef(W_space, xp.GetBlock(1), 0)

    deform = InitialDeformation(dim)
    refconfig = ReferenceConfiguration(dim)

    x_gf.ProjectCoefficient(deform)
    x_ref.ProjectCoefficient(refconfig)
    p_gf.Assign(0.0)

    #  10. Initialize the incompressible neo-Hookean operator
    oper = RubberOperator(spaces, ess_bdr, block_offsets,
                          newton_rel_tol, newton_abs_tol, newton_iter,
                          mu)
    #  11. Solve the Newton system
    try:
        oper.Solve(xp)
        solution = xp.GetDataArray().copy()
    finally:
        oper.shutdown()
        deform.close()
        refconfig.close()

    #  12. Compute the final deformation
    mfem.Subtract(x_gf, x_ref, x_def)

    #  13. Visualize the results if requested
    if (visualization):
        vis_u = mfem.socketstream("localhost", 19916)
        visualize(vis_u, mesh, x_gf, x_def, "Deformation", True)
        vis_p = mfem.socketstream("localhost", 19916)
        visualize(vis_p, mesh, x_gf, p_gf, "Deformation", True)

    #  14. Save the displaced mesh, the final deformation, and the pressure
    nodes = x_gf
    owns_nodes = 0
    nodes, owns_nodes = mesh.SwapNodes(nodes, owns_nodes)

    mesh.Print('deformed.mesh', 8)
    p_gf.Save('pressure.sol', 8)
    x_def.Save("deformation.sol",  8)
    return solution

def visualize(out, mesh, deformed_nodes, field, field_name, init_vis):
    if out is None:
        return

    nodes = deformed_nodes
    owns_nodes = 0

    nodes, owns_nodes = mesh.SwapNodes(nodes, owns_nodes)

    out.precision(8)
    out.send_solution(mesh, field)

    nodes, owns_nodes = mesh.SwapNodes(nodes, owns_nodes)

    if (init_vis):
        out.send_text("window_size 800 800")
        out.send_text("window_title '" + field_name)
        if (mesh.SpaceDimension() == 2):
            out.send_text("view 0 0")
            out.send_text("keys jlA")
        out.send_text("keys cm")         # show colorbar and mesh
        # update value-range; keep mesh-extents fixed
        out.send_text("autoscale value")
    out.flush()


def run(meshfile="beam-quad.mesh", ref_levels=0, order=2, rel_tol=1e-4,
        abs_tol=1e-6, max_iter=500, mu=1.0, visualization=False):
    return ex19_main(SimpleNamespace(mesh=str(meshfile), refine=ref_levels,
                    order=order, relative_tolerance=rel_tol,
                    absolute_tolerance=abs_tol, newton_iterations=max_iter,
                    shear_modulus=mu, visualization=visualization))


if __name__ == "__main__":
    parser = ArgParser(description='Ex19')
    parser.add_argument('-m', '--mesh',
                        default='beam-tet.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-r', '--refine',
                        action='store', default=0, type=int,
                        help="Number of times to refine the mesh uniformly")
    parser.add_argument('-o', '--order',
                        action='store', default=2, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-rel', '--relative-tolerance',
                        action='store', default=1e-4, type=float,
                        help="Relative tolerance for the Newton solve.")
    parser.add_argument('-abs', '--absolute-tolerance',
                        action='store', default=1e-6, type=float,
                        help="Absolute tolerance for the Newton solve.")
    parser.add_argument("-it", "--newton-iterations",
                        action='store', default=500, type=int,
                        help="Maximum iterations for the Newton solve.")
    parser.add_argument("-mu", "--shear-modulus",
                        action='store', default=1.0, type=float,
                        help="Shear modulus in the Neo-Hookean material.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true', default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)
    ex19_main(args)
