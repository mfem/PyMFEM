'''
   multidomain.py

   See c++ version in the MFEM library for more detail
'''

from mfem.par import intArray, doubleArray
import mfem.par as mfem
import os
import sys
from os.path import expanduser, join

from numpy import exp, sqrt, zeros, abs

try:
    from numba import jit
except ImportError:
    assert False, "This example requires numba to run"

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


class ConvectionDiffusionTDO(mfem.PyTimeDependentOperator):
    def __init__(self, fes, ess_tdofs, alpha=1.0, kappa=1e-1):

        mfem.PyTimeDependentOperator.__init__(self, fes.GetTrueVSize())

        # Mass form
        self.Mform = mfem.ParBilinearForm(fes)
        # Stiffness form. Might include diffusion, convection or both.
        self.Kform = mfem.ParBilinearForm(fes)
        # RHS form
        self.bform = mfem.ParLinearForm(fes)
        # Essential true dof array. Relevant for eliminating boundary conditions
        # when using an H1 space.
        self.ess_tdofs_ = intArray(ess_tdofs)
        # Mass matrix solver
        self.M_solver = mfem.CGSolver(fes.GetComm())

        # Mass opeperator
        self.M = mfem.OperatorHandle()
        # Stiffness opeperator. Might include diffusion, convection or both.
        self.K = mfem.OperatorHandle()

        d = mfem.ConstantCoefficient(-kappa)

        @mfem.jit.vector(vdim=fes.GetParMesh().Dimension())
        def velocity_profile(c):
            q = zeros(3)
            A = 1.0
            x = c[0]
            y = c[1]
            r = sqrt(x**2 + y**2)

            if abs(r) >= (0.25 - 1e-8):
                q[2] = 0.0
            else:
                q[2] = A * exp(-x**2 / 2.0 + y**2 / 2.0)
            return q

        q = velocity_profile

        self.Mform.AddDomainIntegrator(mfem.MassIntegrator())
        self.Mform.Assemble(0)
        self.Mform.Finalize()

        if fes.IsDGSpace():
            self.M.Reset(self.Mform.ParallelAssemble(), True)
            inflow = mfem.ConstantCoefficient(0.0)
            self.bform.AddBdrFaceIntegrator(
                mfem.BoundaryFlowIntegrator(inflow, q, alpha))
        else:
            self.Kform.AddDomainIntegrator(
                mfem.ConvectionIntegrator(q, -alpha))
            self.Kform.AddDomainIntegrator(mfem.DiffusionIntegrator(d))
            self.Kform.Assemble(0)

            empty = intArray()
            self.Kform.FormSystemMatrix(empty, self.K)
            self.Mform.FormSystemMatrix(self.ess_tdofs_, self.M)

            self.bform.Assemble()
            self.b = self.bform.ParallelAssemble()

        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(1e-8)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(100)
        self.M_solver.SetPrintLevel(0)

        self.M_prec = mfem.HypreSmoother()
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)

        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.M.Ptr())

        # note we use Height instead of height. SWIG does not
        # allow accessing protected member from Python subclass.
        self.t1 = mfem.Vector(self.Height())
        self.t2 = mfem.Vector(self.Height())

    def Mult(self, u, du_dt):
        self.K.Mult(u, self.t1)
        self.t1 += self.b
        self.M_solver.Mult(self.t1, du_dt)
        du_dt.SetSubVector(self.ess_tdofs_, 0.0)


def run(dt=1e-5,
        t_final=5,
        order=2,
        visualization=True,
        vis_steps=10):

    serial_mesh = mfem.Mesh("multidomain-hex.mesh")
    parent_mesh = mfem.ParMesh(MPI.COMM_WORLD, serial_mesh)

    del serial_mesh

    parent_mesh.UniformRefinement()

    fec = mfem.H1_FECollection(order, parent_mesh.Dimension())

    # Create the sub-domains and accompanying Finite Element spaces from
    # corresponding attributes. This specific mesh has two domain attributes and
    # 9 boundary attributes.
    cylinder_domain_attributes = intArray([1])

    cylinder_submesh = mfem.ParSubMesh.CreateFromDomain(parent_mesh,
                                                        cylinder_domain_attributes)

    print(cylinder_submesh)
    fes_cylinder = mfem.ParFiniteElementSpace(cylinder_submesh, fec)

    inflow_attributes = intArray([0] *
                                 cylinder_submesh.bdr_attributes.Max())
    inflow_attributes[7] = 1

    inner_cylinder_wall_attributes = intArray([0] *
                                              cylinder_submesh.bdr_attributes.Max())
    inner_cylinder_wall_attributes[8] = 1

    # For the convection-diffusion equation inside the cylinder domain, the
    # inflow surface and outer wall are treated as Dirichlet boundary
    # conditions.
    inflow_tdofs = intArray()
    interface_tdofs = intArray()
    ess_tdofs = intArray()
    fes_cylinder.GetEssentialTrueDofs(inflow_attributes,
                                      inflow_tdofs)
    fes_cylinder.GetEssentialTrueDofs(inner_cylinder_wall_attributes,
                                      interface_tdofs)
    ess_tdofs.Append(inflow_tdofs)
    ess_tdofs.Append(interface_tdofs)
    ess_tdofs.Sort()
    ess_tdofs.Unique()

    cd_tdo = ConvectionDiffusionTDO(fes_cylinder, ess_tdofs)

    temperature_cylinder_gf = mfem.ParGridFunction(fes_cylinder)
    temperature_cylinder_gf.Assign(0.0)

    temperature_cylinder = mfem.Vector()
    temperature_cylinder_gf.GetTrueDofs(temperature_cylinder)

    cd_ode_solver = mfem.RK3SSPSolver()
    cd_ode_solver.Init(cd_tdo)

    outer_domain_attributes = intArray([2])

    block_submesh = mfem.ParSubMesh.CreateFromDomain(parent_mesh,
                                                     outer_domain_attributes)

    fes_block = mfem.ParFiniteElementSpace(block_submesh, fec)

    block_wall_attributes = mfem.intArray(
        [0]*block_submesh.bdr_attributes.Max())
    block_wall_attributes[0] = 1
    block_wall_attributes[1] = 1
    block_wall_attributes[2] = 1
    block_wall_attributes[3] = 1

    outer_cylinder_wall_attributes = intArray([0] *
                                              block_submesh.bdr_attributes.Max())
    outer_cylinder_wall_attributes[8] = 1

    fes_block.GetEssentialTrueDofs(block_wall_attributes, ess_tdofs)

    d_tdo = ConvectionDiffusionTDO(fes_block, ess_tdofs, 0.0, 1.0)

    temperature_block_gf = mfem.ParGridFunction(fes_block)
    temperature_block_gf.Assign(0.0)

    one = mfem.ConstantCoefficient(1.0)
    temperature_block_gf.ProjectBdrCoefficient(one, block_wall_attributes)

    temperature_block = mfem.Vector()
    temperature_block_gf.GetTrueDofs(temperature_block)

    d_ode_solver = mfem.RK3SSPSolver()
    d_ode_solver.Init(d_tdo)

    cylinder_surface_attributes = intArray([9])

    cylinder_surface_submesh = mfem.ParSubMesh.CreateFromBoundary(parent_mesh,
                                                                  cylinder_surface_attributes)

    if visualization:
        cyl_sol_sock = mfem.socketstream("localhost", 19916)
        cyl_sol_sock.precision(8)
        cyl_sol_sock << "parallel " << num_procs << " " << myid << "\n"
        cyl_sol_sock.send_solution(cylinder_submesh, temperature_cylinder_gf)

        block_sol_sock = mfem.socketstream("localhost", 19916)
        block_sol_sock.precision(8)
        block_sol_sock << "parallel " << num_procs << " " << myid << "\n"
        block_sol_sock.send_solution(block_submesh, temperature_block_gf)

    temperature_block_to_cylinder_map = mfem.ParSubMesh.CreateTransferMap(
        temperature_block_gf,
        temperature_cylinder_gf)

    t = 0.0
    ti = 1
    last_step = False
    while not last_step:
        if t + dt >= t_final - dt/2:
            last_step = True

        # Advance the diffusion equation on the outer block to the next time step
        t, dt = d_ode_solver.Step(temperature_block, t, dt)

        # Transfer the solution from the inner surface of the outer block to
        # the cylinder outer surface to act as a boundary condition.
        temperature_block_gf.SetFromTrueDofs(temperature_block)

        temperature_block_to_cylinder_map.Transfer(temperature_block_gf,
                                                   temperature_cylinder_gf)

        temperature_cylinder_gf.GetTrueDofs(temperature_cylinder)

        # Advance the convection-diffusion equation on the outer block to the
        # next time step
        t, dt = cd_ode_solver.Step(temperature_cylinder, t, dt)

        if last_step or ti % vis_steps == 0:
            if myid == 0:
                print("step " + str(ti) + ", t = " + str(t))

            temperature_cylinder_gf.SetFromTrueDofs(temperature_cylinder)
            temperature_block_gf.SetFromTrueDofs(temperature_block)

            if visualization:
                cyl_sol_sock << "parallel " << num_procs << " " << myid << "\n"
                cyl_sol_sock << "solution\n" << cylinder_submesh << temperature_cylinder_gf
                cyl_sol_sock.flush()

                block_sol_sock << "parallel " << num_procs << " " << myid << "\n"
                block_sol_sock << "solution\n" << block_submesh << temperature_block_gf
                block_sol_sock.flush()
        ti += 1


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Multidomain')

    parser.add_argument('-o', '--order',
                        action='store', default=2, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-tf', '--t-final',
                        action='store', default=5.0, type=float,
                        help="Final time; start time is 0.")
    parser.add_argument('-dt', '--time-step',
                        action='store', default=1e-5, type=float,
                        help="Time step")
    parser.add_argument('-vis', '--visualization',
                        action='store_true', default=True,
                        help='Enable GLVis visualization')
    parser.add_argument("-vs", "--visualization-steps",
                        action='store', default=10,  type=int,
                        help="Visualize every n-th timestep.")

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    order = args.order
    visualization = args.visualization
    t_final = args.t_final
    dt = args.time_step

    run(dt=dt,
        t_final=t_final,
        order=order,
        visualization=visualization,
        vis_steps=args.visualization_steps)
