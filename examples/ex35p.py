'''
   MFEM example 35p
      See c++ version in the MFEM library for more detail 

   Sample runs:  mpirun -np 4 python ex35p.py -p 0 -o 2
                 mpirun -np 4 python ex35p.py -p 0 -o 2 -pbc '22 23 24' -em 0
                 mpirun -np 4 python ex35p.py -p 1 -o 1 -rp 2
                 mpirun -np 4 python ex35p.py -p 1 -o 2
                 mpirun -np 4 python ex35p.py -p 2 -o 1 -rp 2 -c 15
'''
import os
import sys
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)

mixed = True
visualization = True
pa = False
slu_solver = False  # SuperLU option not supported


def run(mesh_file="",
        order=1,
        ser_ref_levels=0,
        par_ref_levels=0,
        prob=0,
        mode=1,
        herm_conv=True,
        epsilon=1.0,
        sigma=2.,
        mu=1.0,
        device='cpu',
        port_bc_attr=None,
        freq=-1.0,):

    mesh_file = "../data/fichera-mixed.mesh"
    if not mixed or pa:
        mesh_file = "../data/fichera.mesh"

    omega = 2*pi
    if freq > 0.0:
        omega = 2.0 * pi * freq

    if (len(port_bc_attr) == 0 and
        (mesh_file == "../data/fichera-mixed.mesh" or
         mesh_file == "../data/fichera.mesh")):
        port_bc_attr = mfem.intArray([7, 8, 11, 12])

    port_bc_attr.Print()
    conv = mfem.ComplexOperator.HERMITIAN if herm_conv else mfem.ComplexOperator.BLOCK_SYMMETRIC
    # 3. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device)
    if myid == 0:
        device.Print()

    # 4. Read the mesh from the given mesh file. We can handle triangular,
    #    quadrilateral, tetrahedral, hexahedral, surface and volume meshes
    #    with the same code.
    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()

    # 5. Refine the serial mesh on all processors to increase the resolution.
    for i in range(ser_ref_levels):
        mesh.UniformRefinement()

    # 6a. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution. Once the
    #    parallel mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    for i in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 6b. Extract a submesh covering a portion of the boundary
    pmesh_port = mfem.ParSubMesh.CreateFromBoundary(pmesh, port_bc_attr)

    # 7a. Define a parallel finite element space on the parallel mesh. Here we
    #     use continuous Lagrange, Nedelec, or Raviart-Thomas finite elements
    #     of the specified order.
    if dim == 1 and prob != 0:
        if myid == 0:
            print(
                "Switching to problem type 0, H1 basis functions, for 1 dimensional mesh.")
        prob = 0

    if prob == 0:
        fec = mfem.H1_FECollection(order, dim)
    elif prob == 1:
        fec = mfem.ND_FECollection(order, dim)
    elif prob == 2:
        fec = mfem.RT_FECollection(order-1, dim)
    else:
        assert False, "unknown problem"
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    size = fespace.GlobalTrueVSize()
    if myid == 0:
        print("Number of finite element unknowns: " + str(size))

    # 7b. Define a parallel finite element space on the sub-mesh. Here we
    #     use continuous Lagrange, Nedelec, or L2 finite elements of
    #     the specified order.
    if prob == 0:
        fec_port = mfem.H1_FECollection(order, dim-1)
    elif prob == 1:
        if dim == 3:
            fec_port = mfem.ND_FECollection(order, dim-1)
        else:
            fec_port = mfem.L2_FECollection(order - 1, dim-1,
                                            mfem.BasisType.GaussLegendre,
                                            mfem.FiniteElement.INTEGRAL)
    elif prob == 2:
        fec_port = mfem.L2_FECollection(order - 1, dim-1,
                                        mfem.BasisType.GaussLegendre,
                                        mfem.FiniteElement.INTEGRAL)
    else:
        assert False, "should not reach here"

    fespace_port = mfem.ParFiniteElementSpace(pmesh_port, fec_port)
    size_port = fespace_port.GlobalTrueVSize()
    if myid == 0:
        print("Number of finite element port BC unknowns: " + str(size_port))

    # 8a. Define a parallel grid function on the SubMesh which will contain
    #     the field to be applied as a port boundary condition.
    port_bc = mfem.ParGridFunction(fespace_port)
    port_bc.Assign(0.0)

    SetPortBC(prob, dim, mode, port_bc)

    # 8b. Save the SubMesh and associated port boundary condition in parallel.
    #     This output can be viewed later using GLVis:
    #     "glvis -np <np> -m port_mesh -g port_mode"
    pmesh_port.Print("port_mesh"+smyid, 8)
    port_bc.Save("port_mode"+smyid, 8)

    # 8c. Send the port bc, computed on the SubMesh, to a GLVis server.
    if visualization and dim == 3:
        port_sock = mfem.socketstream("localhost", 19916)
        port_sock << "parallel " << num_procs << " " << myid << "\n"
        port_sock.precision(8)
        port_sock << "solution\n" << pmesh_port << port_bc
        port_sock << "window_title 'Port BC'"
        port_sock << "window_geometry 0 0 400 350"
        port_sock.flush()

    # 9. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, the boundary conditions are defined
    #    using an eigenmode of the appropriate type computed on the SubMesh.
    ess_tdof_list = mfem.intArray()
    if pmesh.bdr_attributes.Size():
        ess_bdr = mfem.intArray([1]*pmesh.bdr_attributes.Max())
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 10. Set up the parallel linear form b(.) which corresponds to the
    #     right-hand side of the FEM linear system.
    b = mfem.ParComplexLinearForm(fespace, conv)
    b.Assign(0.0)

    # 11a. Define the solution vector u as a parallel complex finite element
    #      grid function corresponding to fespace. Initialize u to equal zero.
    u = mfem.ParComplexGridFunction(fespace)
    u.Assign(0.0)
    pmesh_port.Transfer(port_bc, u.real())

    # 11b. Send the transferred port bc field to a GLVis server.

    full_bc = mfem.ParGridFunction(fespace)
    port_to_full = mfem.ParTransferMap(port_bc, full_bc)

    full_bc.Assign(0.0)
    port_to_full.Transfer(port_bc, full_bc)

    if visualization:
        full_sock = mfem.socketstream("localhost", 19916)
        full_sock << "parallel " << num_procs << " " << myid << "\n"
        full_sock.precision(8)
        full_sock << "solution\n" << pmesh << full_bc
        full_sock << "window_title 'Transferred BC'"
        full_sock << "window_geometry 400 0 400 350"
        full_sock.flush()

    # 12. Set up the parallel sesquilinear form a(.,.) on the finite element
    #     space corresponding to the damped harmonic oscillator operator of the
    #     appropriate type:
    #
    #     0) A scalar H1 field
    #        -Div(a Grad) - omega^2 b + i omega c
    #
    #     1) A vector H(Curl) field
    #        Curl(a Curl) - omega^2 b + i omega c
    #
    #     2) A vector H(Div) field
    #        -Grad(a Div) - omega^2 b + i omega c
    #
    stiffnessCoef = mfem.ConstantCoefficient(1.0/mu)
    massCoef = mfem.ConstantCoefficient(-omega * omega * epsilon)
    lossCoef = mfem.ConstantCoefficient(omega * sigma)
    negMassCoef = mfem.ConstantCoefficient(omega * omega * epsilon)

    a = mfem.ParSesquilinearForm(fespace, conv)

    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)
    if prob == 0:
        a.AddDomainIntegrator(mfem.DiffusionIntegrator(stiffnessCoef),
                              None)
        a.AddDomainIntegrator(mfem.MassIntegrator(massCoef),
                              mfem.MassIntegrator(lossCoef))
    elif prob == 1:
        a.AddDomainIntegrator(mfem.CurlCurlIntegrator(stiffnessCoef),
                              None)
        a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(massCoef),
                              mfem.VectorFEMassIntegrator(lossCoef))
    elif prob == 2:
        a.AddDomainIntegrator(mfem.DivDivIntegrator(stiffnessCoef),
                              None)
        a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(massCoef),
                              mfem.VectorFEMassIntegrator(lossCoef))
    else:
        assert False, "Unknown probm"

    # 13. Assemble the parallel bilinear form and the corresponding linear
    #     system, applying any necessary transformations such as: parallel
    #     assembly, eliminating boundary conditions, applying conforming
    #     constraints for non-conforming AMR, etc.
    a.Assemble()
    A = mfem.OperatorHandle()
    B = mfem.Vector()
    U = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, u, b, A, U, B)

    if myid == 0:
        print("Size of linear system: " + str(2*size))

    if not slu_solver:
        # 14a. Set up the parallel bilinear form for the preconditioner
        #      corresponding to the appropriate operator
        #
        #      0) A scalar H1 field
        #         -Div(a Grad) - omega^2 b + i omega c
        #
        #      1) A vector H(Curl) field
        #         Curl(a Curl) + omega^2 b + i omega c
        #
        #      2) A vector H(Div) field
        #         -Grad(a Div) - omega^2 b + i omega c
        pcOp = mfem.ParBilinearForm(fespace)
        if pa:
            pcOp.SetAssemblyLevel(mfem.AssemblyLevel_PARTIAL)

        if prob == 0:
            pcOp.AddDomainIntegrator(mfem.DiffusionIntegrator(stiffnessCoef))
            pcOp.AddDomainIntegrator(mfem.MassIntegrator(massCoef))
            pcOp.AddDomainIntegrator(mfem.MassIntegrator(lossCoef))
        elif prob == 1:
            pcOp.AddDomainIntegrator(mfem.CurlCurlIntegrator(stiffnessCoef))
            pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(negMassCoef))
            pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(lossCoef))
        elif prob == 2:
            pcOp.AddDomainIntegrator(mfem.DivDivIntegrator(stiffnessCoef))
            pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(massCoef))
            pcOp.AddDomainIntegrator(mfem.VectorFEMassIntegrator(lossCoef))
        pcOp.Assemble()

        # 14b. Define and apply a parallel FGMRES solver for AU=B with a block
        #      diagonal preconditioner based on the appropriate multigrid
        #      preconditioner from hypre.
        blockTrueOffsets = mfem.intArray()
        blockTrueOffsets.SetSize(3)
        blockTrueOffsets[0] = 0
        blockTrueOffsets[1] = A.Height() // 2
        blockTrueOffsets[2] = A.Height() // 2
        blockTrueOffsets.PartialSum()

        BDP = mfem.BlockDiagonalPreconditioner(blockTrueOffsets)

        if pa:
            pc_r = mfem.OperatorJacobiSmoother(pcOp, ess_tdof_list)
            pc_i = None
        else:
            PCOp = mfem.OperatorHandle()
            pcOp.FormSystemMatrix(ess_tdof_list, PCOp)
            if prob == 0:
                pc_r = mfem.HypreBoomerAMG(PCOp.AsHypreParMatrix())
            elif prob == 1:
                pc_r = mfem.HypreAMS(PCOp.AsHypreParMatrix(), fespace)
            elif prob == 2:
                if dim == 2:
                    pc_r = mfem.HypreAMS(PCOp.AsHypreParMatrix(), fespace)
                else:
                    pc_r = mfem.HypreADS(PCOp.AsHypreParMatrix(), fespace)

        pc_i = mfem.ScaledOperator(pc_r,
                                   -1 if conv == mfem.ComplexOperator.HERMITIAN else 1)
        BDP.SetDiagonalBlock(0, pc_r)
        BDP.SetDiagonalBlock(1, pc_i)

        fgmres = mfem.FGMRESSolver(MPI.COMM_WORLD)
        fgmres.SetPreconditioner(BDP)
        fgmres.SetOperator(A.Ptr())
        fgmres.SetRelTol(1e-6)
        fgmres.SetMaxIter(1000)
        fgmres.SetPrintLevel(1)
        fgmres.Mult(B, U)

    # 15. Recover the parallel grid function corresponding to U. This is the
    #     local finite element solution on each processor.
    a.RecoverFEMSolution(U, b, u)

    # 16. Save the refined mesh and the solution in parallel. This output can be
    #     viewed later using GLVis: "glvis -np <np> -m mesh -g sol_r" or
    #     "glvis -np <np> -m mesh -g sol_i".
    pmesh.Print("mesh"+smyid, 8)
    u.real().Save("sol_r"+smyid, 8)
    u.imag().Save("sol_i"+smyid, 8)

    # 17. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock_r = mfem.socketstream("localhost", 19916)
        sol_sock_r << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_r.precision(8)
        sol_sock_r << "solution\n" << pmesh << u.real(
        ) << "window_title 'Solution: Real Part'"
        sol_sock_r << "window_geometry 800 0 400 350"
        sol_sock_r.flush()

        sol_sock_i = mfem.socketstream("localhost", 19916)
        sol_sock_i << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_i.precision(8)
        sol_sock_i << "solution\n" << pmesh << u.imag(
        ) << "window_title 'Solution: Imag Part'"
        sol_sock_i << "window_geometry 1200 0 400 350"
        sol_sock_i.flush()

    if visualization:
        u_t = mfem.ParGridFunction(fespace)
        u_t.Assign(u.real())
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << pmesh << u_t
        sol_sock << "window_title 'Harmonic Solution (t = 0.0 T)'"
        sol_sock << "window_geometry 0 432 600 450"
        sol_sock << "pause\n"
        sol_sock.flush()

        if myid == 0:
            print(
                "GLVis visualization paused. Press space (in the GLVis window) to resume it.")
        num_frames = 32
        i = 0

        # Let's plot one wave cycle...
        for i in range(160):
            t = (i % num_frames) / num_frames
            oss = "Harmonic Solution (t = " + str(t) + " T)"
            dd = (cos(2.0 * pi * t)*u.real().GetDataArray() +
                  sin(2.0 * pi * t)*u.imag().GetDataArray())
            # we can not load numpy directly...(sorry)
            u_t.Assign(mfem.Vector(dd))
            sol_sock << "parallel " << num_procs << " " << myid << "\n"
            sol_sock << "solution\n" << pmesh << u_t
            sol_sock << "window_title '" << oss << "'"
            sol_sock.flush()


'''
   Solves the eigenvalue problem -Div(Grad x) = lambda x with homogeneous
   Dirichlet boundary conditions on the boundary of the domain. Returns mode
   number "mode" (counting from zero) in the ParGridFunction "x".
'''


def ScalarWaveGuide(mode, x):
    nev = max(mode + 2, 5)
    seed = 75

    fespace = x.ParFESpace()
    pmesh = fespace.GetParMesh()

    if pmesh.bdr_attributes.Size() > 0:
        ess_bdr = mfem.intArray([1]*pmesh.bdr_attributes.Max())
    else:
        ess_bdr = mfem.intArray()

    a = mfem.ParBilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator())
    a.Assemble()
    a.EliminateEssentialBCDiag(ess_bdr, 1.0)
    a.Finalize()

    m = mfem.ParBilinearForm(fespace)
    m.AddDomainIntegrator(mfem.MassIntegrator())
    m.Assemble()

    # shift the eigenvalue corresponding to eliminated dofs to a large value
    m.EliminateEssentialBCDiag(ess_bdr, sys.float_info.min)
    m.Finalize()

    A = a.ParallelAssemble()
    M = m.ParallelAssemble()

    amg = mfem.HypreBoomerAMG(A)
    amg.SetPrintLevel(0)

    lobpcg = mfem.HypreLOBPCG(MPI.COMM_WORLD)
    lobpcg.SetNumModes(nev)
    lobpcg.SetRandomSeed(seed)
    lobpcg.SetPreconditioner(amg)
    lobpcg.SetMaxIter(200)
    lobpcg.SetTol(1e-8)
    lobpcg.SetPrecondUsageMode(1)
    lobpcg.SetPrintLevel(1)
    lobpcg.SetMassMatrix(M)
    lobpcg.SetOperator(A)
    lobpcg.Solve()

    x.Assign(lobpcg.GetEigenvector(mode))


'''
   Solves the eigenvalue problem -Curl(Curl x) = lambda x with homogeneous
   Dirichlet boundary conditions, on the tangential component of x, on the
   boundary of the domain. Returns mode number "mode" (counting from zero) in
   the ParGridFunction "x".
'''


def VectorWaveGuide(mode, x):
    nev = max(mode + 2, 5)
    fespace = x.ParFESpace()
    pmesh = fespace.GetParMesh()

    if pmesh.bdr_attributes.Size() > 0:
        ess_bdr = mfem.intArray([1]*pmesh.bdr_attributes.Max())
    else:
        ess_bdr = mfem.intArray()

    a = mfem.ParBilinearForm(fespace)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator())
    a.Assemble()
    a.EliminateEssentialBCDiag(ess_bdr, 1.0)
    a.Finalize()

    m = mfem.ParBilinearForm(fespace)
    m.AddDomainIntegrator(mfem.VectorFEMassIntegrator())
    m.Assemble()

    # shift the eigenvalue corresponding to eliminated dofs to a large value
    m.EliminateEssentialBCDiag(ess_bdr, sys.float_info.min)
    m.Finalize()

    A = a.ParallelAssemble()
    M = m.ParallelAssemble()

    ams = mfem.HypreAMS(A, fespace)
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
    ame.Solve()

    x.Assign(ame.GetEigenvector(mode))


'''
   Solves the eigenvalue problem -Div(Grad x) = lambda x with homogeneous
   Neumann boundary conditions on the boundary of the domain. Returns mode
   number "mode" (counting from zero) in the ParGridFunction "x_l2". Note that
   mode 0 is a constant field so higher mode numbers are often more
   interesting. The eigenmode is solved using continuous H1 basis of the
   appropriate order and then projected onto the L2 basis and returned.
'''


def PseudoScalarWaveGuide(mode, x_l2):
    nev = max(mode + 2, 5)
    seed = 75

    fespace_l2 = x_l2.ParFESpace()
    pmesh = fespace_l2.GetParMesh()
    order_l2 = fespace_l2.FEColl().GetOrder()

    fec = mfem.H1_FECollection(order_l2+1, pmesh.Dimension())
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    x = mfem.ParGridFunction(fespace)
    x.Assign(0.0)

    xCoef = mfem.GridFunctionCoefficient(x)
    if mode == 0:
        x.Assign(1.0)
        x_l2.ProjectCoefficient(xCoef)

    a = mfem.ParBilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator())
    a.AddDomainIntegrator(mfem.MassIntegrator())
    a.Assemble()
    a.Finalize()

    m = mfem.ParBilinearForm(fespace)
    m.AddDomainIntegrator(mfem.MassIntegrator())
    m.Assemble()
    m.Finalize()

    A = a.ParallelAssemble()
    M = m.ParallelAssemble()

    amg = mfem.HypreBoomerAMG(A)
    amg.SetPrintLevel(0)

    lobpcg = mfem.HypreLOBPCG(MPI.COMM_WORLD)
    lobpcg.SetNumModes(nev)
    lobpcg.SetRandomSeed(seed)
    lobpcg.SetPreconditioner(amg)
    lobpcg.SetMaxIter(200)
    lobpcg.SetTol(1e-8)
    lobpcg.SetPrecondUsageMode(1)
    lobpcg.SetPrintLevel(1)
    lobpcg.SetMassMatrix(M)
    lobpcg.SetOperator(A)
    lobpcg.Solve()

    x.Assign(lobpcg.GetEigenvector(mode))
    x_l2.ProjectCoefficient(xCoef)

# Compute eigenmode "mode" of either a Dirichlet or Neumann Laplacian or of a
# Dirichlet curl curl operator based on the problem type and dimension of the
# domain.


def SetPortBC(prob, dim,  mode, port_bc):
    if prob == 0:
        ScalarWaveGuide(mode, port_bc)
    elif prob == 1:
        if dim == 3:
            VectorWaveGuide(mode, port_bc)
        else:
            PseudoScalarWaveGuide(mode, port_bc)
    elif prob == 2:
        PseudoScalarWaveGuide(mode, port_bc)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex35p (Port Boundary Conditions using SubMesh Transfers)')
    parser.add_argument('-m', '--mesh',
                        default='inline-quad.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument("-rs", "--refine-serial",
                        action='store', type=int, default=1,
                        help="Number of times to refine the mesh uniformly in serial.")
    parser.add_argument("-rp", "--refine-parallel",
                        action='store', type=int, default=1,
                        help="Number of times to refine the mesh uniformly in paralle.")
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-p", "--problem-type",
                        action='store', type=int, default=0,
                        help="\n".join(["Choose between 0: H_1, 1: H(Curl), or 2: H(Div) "
                                        "damped harmonic oscillator."]))
    parser.add_argument("-em", "--eigenmode",
                        action='store', type=int, default=1,
                        help="Choose the index of the port eigenmode.")
    parser.add_argument("-a", "--stiffness-coef",
                        action='store', type=float, default=0.0,
                        help="Stiffness coefficient (spring constant or 1/mu).")
    parser.add_argument("-b", "--mass-coef",
                        action='store', type=float, default=1.0,
                        help="Mass coefficient (or epsilon).")
    parser.add_argument("-c", "--damping-coef",
                        action='store', type=float, default=2.0,
                        help="Damping coefficient (or sigma).")
    parser.add_argument("-mu", "--permeability",
                        action='store', type=float, default=1.0,
                        help="Permeability of free space (or 1/(spring constant)).")
    parser.add_argument("-eps", "--permittivity",
                        action='store', type=float, default=1.0,
                        help="Permittivity of free space (or mass constant).")
    parser.add_argument("-sigma", "--conductivity",
                        action='store', type=float, default=2.0,
                        help="Conductivity (or damping constant).")
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=-1.0,
                        help="Set the frequency for the exact")
    parser.add_argument("-pbc", "--port-bc-attr",
                        action='store', type=str, default="",
                        help="Attributes of port boundary condition")
    parser.add_argument("-herm", "--hermitian",
                        action='store_true',
                        default=True,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument("-no-herm", "--no-hermitian",
                        action='store_true',
                        default=False,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument("-no-vis", "--no-visualization",
                        action='store_true', default=False,
                        help='Enable GLVis visualization')
    parser.add_argument("-hex", "--hex-mesh",
                        action='store_true',  default=False,
                        help="Mixed mesh of hexahedral mesh.")
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")

    args = parser.parse_args()

    herm = False if args.no_hermitian else True
    args.hermitian = herm
    args.no_hermitian = not herm

    if myid == 0:
        parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    port_bc_attr = [int(x) for x in args.port_bc_attr.split(' ') if len(x) > 0]
    globals()["mixed"] = not args.hex_mesh
    globals()["pa"] = args.partial_assembly
    globals()["visualization"] = not args.no_visualization

    if args.stiffness_coef != 0.0:
        mu = 1.0/args.stiffness_coef
    elif args.permeability != 1.0:
        mu = args.permeability
    else:
        mu = 1.0

    if args.permittivity != 1.0:
        epsilon = args.permittivity
    elif args.mass_coef != 1.0:
        epsilon = args.mass_coef
    else:
        epsilon = 1.0

    if args.conductivity != 2.0:
        sigma = args.conductivity
    elif args.damping_coef != 2.0:
        sigma = args.damping_coef
    else:
        sigma = 2.0

    run(mesh_file=meshfile,
        order=args.order,
        ser_ref_levels=args.refine_serial,
        par_ref_levels=args.refine_parallel,
        prob=args.problem_type,
        mode=args.eigenmode,
        herm_conv=herm,
        epsilon=epsilon,
        sigma=sigma,
        mu=mu,
        device=args.device,
        port_bc_attr=port_bc_attr,
        freq=args.frequency)
