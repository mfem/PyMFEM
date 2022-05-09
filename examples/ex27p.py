'''
   MFEM example 27p
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array, floor, log, arcsin

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def generate_serial_mesh(ref, a_):

    mesh = mfem.Mesh(2, 29, 16, 24, 2)

    for i in range(2):
        o = 13 * i
        mesh.AddQuad([o + 0, o + 3, o + 4, o + 1])
        mesh.AddQuad([o + 1, o + 4, o + 5, o + 2])
        mesh.AddQuad([o + 5, o + 8, o + 9, o + 2])
        mesh.AddQuad([o + 8, o + 12, o + 15, o + 9])
        mesh.AddQuad([o + 11, o + 14, o + 15, o + 12])
        mesh.AddQuad([o + 10, o + 13, o + 14, o + 11])
        mesh.AddQuad([o + 6, o + 13, o + 10, o + 7])
        mesh.AddQuad([o + 0, o + 6, o + 7, o + 3])

    mesh.AddBdrSegment([0, 6], 1)
    mesh.AddBdrSegment([6, 13], 1)
    mesh.AddBdrSegment([13, 19], 1)
    mesh.AddBdrSegment([19, 26], 1)

    mesh.AddBdrSegment([28, 22], 2)
    mesh.AddBdrSegment([22, 15], 2)
    mesh.AddBdrSegment([15, 9], 2)
    mesh.AddBdrSegment([9, 2], 2)

    for i in range(2):
        o = 13 * i
        mesh.AddBdrSegment([o+7, o+3], 3 + i)
        mesh.AddBdrSegment([o+10, o+7], 3 + i)
        mesh.AddBdrSegment([o+11, o+10], 3 + i)
        mesh.AddBdrSegment([o+12, o+11], 3 + i)
        mesh.AddBdrSegment([o+8, o+12], 3 + i)
        mesh.AddBdrSegment([o+5, o+8], 3 + i)
        mesh.AddBdrSegment([o+4, o+5], 3 + i)
        mesh.AddBdrSegment([o+3, o+4], 3 + i)

    a = a_ / sqrt(2)
    mesh.AddVertex([-1.0, -0.5])
    mesh.AddVertex([-1.0, 0.0])
    mesh.AddVertex([-1.0, 0.5])

    mesh.AddVertex([-0.5 - a, -a])
    mesh.AddVertex([-0.5 - a, 0.0])
    mesh.AddVertex([-0.5 - a, a])

    mesh.AddVertex([-0.5, -0.5])
    mesh.AddVertex([-0.5, -a])
    mesh.AddVertex([-0.5, a])
    mesh.AddVertex([-0.5, 0.5])

    mesh.AddVertex([-0.5 + a, -a])
    mesh.AddVertex([-0.5 + a, 0.0])
    mesh.AddVertex([-0.5 + a, a])

    mesh.AddVertex([.0, -0.5])
    mesh.AddVertex([.0, 0.0])
    mesh.AddVertex([.0, 0.5])

    mesh.AddVertex([0.5 - a, -a])
    mesh.AddVertex([0.5 - a, 0.0])
    mesh.AddVertex([0.5 - a, a])

    mesh.AddVertex([0.5, -0.5])
    mesh.AddVertex([0.5, -a])
    mesh.AddVertex([0.5, a])
    mesh.AddVertex([0.5, 0.5])

    mesh.AddVertex([0.5 + a, -a])
    mesh.AddVertex([0.5 + a, 0.0])
    mesh.AddVertex([0.5 + a, a])

    mesh.AddVertex([1.0, -0.5])
    mesh.AddVertex([1.0, 0.0])
    mesh.AddVertex([1.0, 0.5])

    mesh.FinalizeTopology()

    mesh.SetCurvature(1, True)

    # Stitch the ends of the stack together
    # In Python, we just fill list like this
    v2v = [i for i in range(mesh.GetNV() - 3)]
    v2v.append(0)
    v2v.append(1)
    v2v.append(2)

    # renumber elements
    for i in range(mesh.GetNE()):
        el = mesh.GetElement(i)
        nv = el.GetNVertices()
        # We need to re-write the vertex data.
        # el.GetVertices returns *int. So we put it into intArray
        # to access its element
        v = mfem.intArray([el.GetVertices(), nv])
        for j in range(nv):
            v[j] = v2v[v[j]]

          #  renumber boundary elements
        for i in range(mesh.GetNBE()):
            el = mesh.GetBdrElement(i)
            nv = el.GetNVertices()
            v = mfem.intArray([el.GetVertices(), nv])
            for j in range(nv):
                v[j] = v2v[v[j]]

    mesh.RemoveUnusedVertices()
    mesh.RemoveInternalBoundaries()

    mesh.SetCurvature(3, True)

    for l in range(ref):
        mesh.UniformRefinement()

    sdim = mesh.SpaceDimension()

    def quad_trans(u, v, write=False):
        a = a_
        d = 4.0 * a * (sqrt(2) - 2 * a) * (1.0 - 2.0 * v)
        v0 = ((1.0 + sqrt(2)) * (sqrt(2) * a - 2.0 * v) *
              ((4.0 - 3 * sqrt(2)) * a +
               (8.0 * (sqrt(2) - 1.0) * a - 2.0) * v) / d)

        r = (2.0 * ((sqrt(2) - 1.0) * a * a * (1.0 - 4.0 * v) +
                    2.0 * (1.0 + sqrt(2) *
                           (1.0 + 2.0 * (2.0 * a - sqrt(2) - 1.0) * a)) * v * v
                    ) / d)

        t = arcsin(v / r) * u / v
        if write:
            print("u, v, r, v0, t " +
                  "{:g}".format(u) + " " +
                  "{:g}".format(v) + " " +
                  "{:g}".format(r) + " " +
                  "{:g}".format(v0) + " " +
                  "{:g}".format(t))

        x = r * sin(t)
        y = r * cos(t) - v0

        return x, y

    class cTrans(mfem.VectorPyCoefficient):
        def __init__(self):
            mfem.VectorPyCoefficient.__init__(self, sdim)

        def EvalValue(self, u):
            tol = 1e-4
            x = u*0
            if u[1] > 0.5 - tol or u[1] < -0.5 + tol:
                return u
            if (u[0] > 1.0 - tol or u[0] < -1.0 + tol or abs(u[0]) < tol):
                return u
            if u[0] > 0.0:
                if u[1] > abs(u[0] - 0.5):
                    x0, x1 = quad_trans(u[0] - 0.5, u[1])
                    x[0] = x0 + 0.5
                    x[1] = x1
                    return x
                if u[1] < -abs(u[0] - 0.5):
                    x0, x1 = quad_trans(u[0] - 0.5, -u[1])
                    x[0] = x0 + 0.5
                    x[1] = -x1
                    return x
                if u[0] - 0.5 > abs(u[1]):
                    x1, x0 = quad_trans(u[1], u[0] - 0.5)
                    x[0] = x0+0.5
                    x[1] = x1
                    return x
                if u[0] - 0.5 < -abs(u[1]):
                    x1, x0 = quad_trans(u[1], 0.5 - u[0])
                    x[0] = -x0
                    x[0] = x[0] + 0.5
                    x[1] = x1
                    return x
                else:
                    pass
            else:
                if u[1] > abs(u[0] + 0.5):
                    x0, x1 = quad_trans(u[0] + 0.5, u[1])
                    x[0] = x0 - 0.5
                    x[1] = x1
                    return x
                if u[1] < -abs(u[0] + 0.5):
                    x0, x1 = quad_trans(u[0] + 0.5, -u[1])
                    x[0] = x0 - 0.5
                    x[1] = -x1
                    return x
                if u[0] + 0.5 > abs(u[1]):
                    x1, x0 = quad_trans(u[1], u[0] + 0.5)
                    x[0] = x0 - 0.5
                    x[1] = x1
                    return x
                if u[0] + 0.5 < -abs(u[1]):
                    x1, x0 = quad_trans(u[1], -0.5 - u[0])
                    x[0] = -x0
                    x[0] = x[0] - 0.5
                    x[1] = x1
                    return x
            x = u
            return x

    mesh.Transform(cTrans())

    return mesh


def IntegrateBC(x, bdr, alpha, beta, gamma):
    nrm = 0.0
    avg = 0.0
    err = 0.0

    a_is_zero = (alpha == 0.0)
    b_is_zero = (beta == 0.0)

    fes = x.ParFESpace()
    assert fes.GetVDim() == 1, ""

    mesh = fes.GetParMesh()
    shape = mfem.Vector()
    loc_dofs = mfem.Vector()
    w_nor = mfem.Vector()
    dshape = mfem.DenseMatrix()

    dof_ids = mfem.intArray()

    battrs = mesh.GetBdrAttributeArray()

    for i in range(mesh.GetNBE()):
        if bdr[battrs[i]-1] == 0:
            continue

        FTr = mesh.GetBdrFaceTransformations(i)
        if FTr is None:
            continue

        fe = fes.GetFE(FTr.Elem1No)
        assert fe.GetMapType() == mfem.FiniteElement.VALUE, ""

        int_order = 2*fe.GetOrder() + 3

        ir = mfem.IntRules.Get(FTr.GetGeometryType(), int_order)

        dof_ids = fes.GetElementDofs(FTr.Elem1No)

        x.GetSubVector(mfem.intArray(dof_ids), loc_dofs)
        if not a_is_zero:
            sdim = FTr.Face.GetSpaceDim()
            w_nor.SetSize(sdim)
            dshape.SetSize(fe.GetDof(), sdim)
        if not b_is_zero:
            shape.SetSize(fe.GetDof())

        for j in range(ir.GetNPoints()):
            ip = ir.IntPoint(j)
            eip = mfem.IntegrationPoint()
            FTr.Loc1.Transform(ip, eip)
            FTr.Face.SetIntPoint(ip)
            face_weight = FTr.Face.Weight()
            val = 0.0
            if not a_is_zero:
                FTr.Elem1.SetIntPoint(eip)
                fe.CalcPhysDShape(FTr.Elem1, dshape)
                mfem.CalcOrtho(FTr.Face.Jacobian(), w_nor)
                val += alpha * \
                    dshape.InnerProduct(w_nor, loc_dofs) / face_weight
            if not b_is_zero:
                fe.CalcShape(eip, shape)
                val += beta * (shape * loc_dofs)

            # Measure the length of the boundary
            nrm += ip.weight * face_weight

            # Integrate alpha * n.Grad(x) + beta * x
            avg += val * ip.weight * face_weight

            # Integrate |alpha * n.Grad(x) + beta * x - gamma|^2
            val -= gamma
            err += (val*val) * ip.weight * face_weight

        # Normalize by the length of the boundary
    glb_nrm = MPI.COMM_WORLD.allreduce(nrm, op=MPI.SUM)
    glb_avg = MPI.COMM_WORLD.allreduce(avg, op=MPI.SUM)
    glb_err = MPI.COMM_WORLD.allreduce(err, op=MPI.SUM)

    if abs(glb_nrm) > 0.0:
        glb_err /= glb_nrm
        glb_avg /= glb_nrm

    # Compute l2 norm of the error in the boundary condition (negative
    # quadrature weights may produce negative 'err')
    glb_err = sqrt(abs(glb_err))

    # Return the average value of alpha * n.Grad(x) + beta * x
    return glb_err, glb_avg


def run(order=1,
        h1=True,
        sigma=-1.0,
        kappa=-1.0,
        ser_ref_levels=2,
        par_ref_levels=1,
        mat_val=1.0,
        dbc_val=0.0,
        nbc_val=1.0,
        rbc_a_val=1.0,
        rbc_b_val=1.0,
        a_=0.2,
        visualization=True):

    device = mfem.Device('cpu')
    if myid == 0: device.Print()

    if kappa < 0 and not h1:
        kappa = (order+1.0)**2
    if a_ < 0.01:
        print("Hole radius too small, resetting to 0.01.")
        a_ = 0.01
    if a_ > 0.49:
        print("Hole radius too large, resetting to 0.49.")
        a_ = 0.49

    # 2. Construct the (serial) mesh and refine it if requested.
    mesh = generate_serial_mesh(ser_ref_levels, a_)
    dim = mesh.Dimension()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    for i in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 3. Define a finite element space on the serial mesh. Here we use either
    #    continuous Lagrange finite elements or discontinuous Galerkin finite
    #   elements of the specified order.
    fec = (mfem.H1_FECollection(order, dim) if h1 else
           mfem.DG_FECollection(order, dim))
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    size = fespace.GlobalTrueVSize()
    if myid == 0:
        print("Number of finite element unknowns: " + str(size))

    # 4. Create "marker arrays" to define the portions of boundary associated
    #    with each type of boundary condition. These arrays have an entry
    #    corresponding to each boundary attribute.  Placing a '1' in entry i
    #    marks attribute i+1 as being active, '0' is inactive.
    nbc_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
    rbc_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
    dbc_bdr = mfem.intArray(pmesh.bdr_attributes.Max())

    nbc_bdr.Assign(0)
    nbc_bdr[0] = 1
    rbc_bdr.Assign(0)
    rbc_bdr[1] = 1
    dbc_bdr.Assign(0)
    dbc_bdr[2] = 1

    ess_tdof_list = mfem.intArray()
    if h1 and pmesh.bdr_attributes.Size():
        # For a continuous basis the linear system must be modified to enforce an
        # essential (Dirichlet) boundary condition. In the DG case this is not
        # necessary as the boundary condition will only be enforced weakly.
        fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list)

    # 5. Setup the various coefficients needed for the Laplace operator and the
    #    various boundary conditions. In general these coefficients could be
    #    functions of position but here we use only constants.
    matCoef = mfem.ConstantCoefficient(mat_val)
    dbcCoef = mfem.ConstantCoefficient(dbc_val)
    nbcCoef = mfem.ConstantCoefficient(nbc_val)
    rbcACoef = mfem.ConstantCoefficient(rbc_a_val)
    rbcBCoef = mfem.ConstantCoefficient(rbc_b_val)

    # Since the n.Grad(u) terms arise by integrating -Div(m Grad(u)) by parts we
    # must introduce the coefficient 'm' into the boundary conditions.
    # Therefore, in the case of the Neumann BC, we actually enforce m n.Grad(u)
    # = m g rather than simply n.Grad(u) = g.
    m_nbcCoef = mfem.ProductCoefficient(matCoef, nbcCoef)
    m_rbcACoef = mfem.ProductCoefficient(matCoef, rbcACoef)
    m_rbcBCoef = mfem.ProductCoefficient(matCoef, rbcBCoef)

    # 6. Define the solution vector u as a finite element grid function
    #    corresponding to fespace. Initialize u with initial guess of zero.
    u = mfem.ParGridFunction(fespace)
    u.Assign(0.0)

    # 7. Set up the bilinear form a(.,.) on the finite element space
    #    corresponding to the Laplacian operator -Delta, by adding the Diffusion
    #    domain integrator.
    a = mfem.ParBilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(matCoef))
    if h1:
        # Add a Mass integrator on the Robin boundary
        a.AddBoundaryIntegrator(mfem.MassIntegrator(m_rbcACoef), rbc_bdr)
    else:
        # Add the interfacial portion of the Laplace operator
        a.AddInteriorFaceIntegrator(mfem.DGDiffusionIntegrator(matCoef,
                                                               sigma, kappa))

        # Counteract the n.Grad(u) term on the Dirichlet portion of the boundary
        a.AddBdrFaceIntegrator(mfem.DGDiffusionIntegrator(matCoef, sigma, kappa),
                               dbc_bdr)

        # Augment the n.Grad(u) term with a*u on the Robin portion of boundary
        a.AddBdrFaceIntegrator(
            mfem.BoundaryMassIntegrator(m_rbcACoef), rbc_bdr)
    a.Assemble()

    # 8. Assemble the linear form for the right hand side vector.
    b = mfem.ParLinearForm(fespace)

    if h1:
        # Set the Dirichlet values in the solution vector
        u.ProjectBdrCoefficient(dbcCoef, dbc_bdr)

        # Add the desired value for n.Grad(u) on the Neumann boundary
        b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(m_nbcCoef), nbc_bdr)

        # Add the desired value for n.Grad(u) + a*u on the Robin boundary
        b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(m_rbcBCoef), rbc_bdr)
    else:
        # Add the desired value for the Dirichlet boundary
        b.AddBdrFaceIntegrator(mfem.DGDirichletLFIntegrator(dbcCoef, matCoef,
                                                            sigma, kappa),
                               dbc_bdr)

        # Add the desired value for n.Grad(u) on the Neumann boundary
        b.AddBdrFaceIntegrator(mfem.BoundaryLFIntegrator(m_nbcCoef), nbc_bdr)

        # Add the desired value for n.Grad(u) + a*u on the Robin boundary
        b.AddBdrFaceIntegrator(mfem.BoundaryLFIntegrator(m_rbcBCoef), rbc_bdr)

    b.Assemble()

    # 9. Construct the linear system.
    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, u, b, A, X, B)

    # 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
    #     solve the system AX=B with PCG in the symmetric case, and GMRES in the
    #     non-symmetric one.
    amg = mfem.HypreBoomerAMG()
    if sigma == -1.0:
        pcg = mfem.HyprePCG(MPI.COMM_WORLD)
        pcg.SetTol(1e-12)
        pcg.SetMaxIter(200)
        pcg.SetPrintLevel(2)
        pcg.SetPreconditioner(amg)
        pcg.SetOperator(A.Ptr())
        pcg.Mult(B, X)
    else:
        gmres.GMRESSolver(MPI.COMM_WORLD)
        gmres.SetAbsTol(0.0)
        gmres.SetRelTol(1e-12)
        gmres.SetMaxIter(200)
        gmres.SetKDim(10)
        gmres.SetPrintLevel(1)
        gmres.SetPreconditioner(amg)
        gmres.SetOperator(A.Ptr())
        gmres.Mult(B, X)

    # 12. Recover the grid function corresponding to U. This is the local finite
    #     element solution.
    a.RecoverFEMSolution(X, b, u)

    # 13. Build a mass matrix to help solve for n.Grad(u) where 'n' is a surface
    #     normal.
    m = mfem.ParBilinearForm(fespace)
    m.AddDomainIntegrator(mfem.MassIntegrator())
    m.Assemble()

    ess_tdof_list.SetSize(0)
    M = mfem.OperatorPtr()
    m.FormSystemMatrix(ess_tdof_list, M)

    # 14. Compute the various boundary integrals.
    if myid == 0:
        print("Verifying boundary conditions" +
              "=============================")

    # Integrate the solution on the Dirichlet boundary and compare to the
    # expected value.
    def print0(*args):
        if myid == 0:
            print(*args)

    err, avg = IntegrateBC(u, dbc_bdr, 0.0, 1.0, dbc_val)

    hom_dbc = (dbc_val == 0.0)
    err /= 1.0 if hom_dbc else abs(dbc_val)
    print0("Average of solution on Gamma_dbc:\t" +
           "{:g}".format(avg) + ", \t" +
           ("absolute" if hom_dbc else "relative") +
           " error " + "{:g}".format(err))

    # Integrate n.Grad(u) on the inhomogeneous Neumann boundary and compare
    # to the expected value.
    err, avg = IntegrateBC(u, nbc_bdr, 1.0, 0.0, nbc_val)

    hom_nbc = (nbc_val == 0.0)
    err /= 1.0 if hom_nbc else abs(nbc_val)
    print0("Average of n.Grad(u) on Gamma_ndbc:\t" +
           "{:g}".format(avg) + ", \t" +
           ("absolute" if hom_nbc else "relative") +
           " error " + "{:g}".format(err))

    # Integrate n.Grad(u) on the homogeneous Neumann boundary and compare to
    # the expected value of zero.
    nbc0_bdr = mfem.intArray(mesh.bdr_attributes.Max())
    nbc0_bdr.Assign(0)
    nbc0_bdr[3] = 1

    err, avg = IntegrateBC(u, nbc0_bdr, 1.0, 0.0, 0.0)
    hom_nbc = True
    print0("Average of n.Grad(u) on Gamma_ndbc0:\t" +
           "{:g}".format(avg) + ", \t" +
           ("absolute" if hom_nbc else "relative") +
           " error " + "{:g}".format(err))

    # Integrate n.Grad(u) + a * u on the Robin boundary and compare to the
    # expected value.
    err, avg = IntegrateBC(u, rbc_bdr, 1.0, rbc_a_val, rbc_b_val)

    hom_rbc = (rbc_b_val == 0.0)
    err /= 1.0 if hom_rbc else abs(nbc_val)
    print0("Average of n.Grad(u)+a*u on Gamma_rdbc:\t" +
           "{:g}".format(avg) + ", \t" +
           ("absolute" if hom_rbc else "relative") +
           " error " + "{:g}".format(err))

    # 15. Save the refined mesh and the solution. This output can be viewed
    #     later using GLVis: "glvis -m refined.mesh -g sol.gf".
    pmesh.Print("mesh"+smyid, 8)
    u.Save("sol"+smyid, 8)

    # 16. Send the solution by socket to a GLVis server.
    if visualization:
        title_str = "H1" if h1 else "DG"
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << u
        sol_sock << "window_title '" << title_str << " Solution'"
        sol_sock << " keys 'mmc'"
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex27 (Laplace boundary conditionss)')

    parser.add_argument("-h1", "--continuous",
                        action='store_true', default=True,
                        help='Select continuous "H1" element')
    parser.add_argument("-dg", "--discontinuous",
                        action='store_true', default=False,
                        help='Select continuous "DG" element')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-s", "--sigma",
                        action='store', default=-1.0, type=float,
                        help="One of the two DG penalty parameters, typically +1/-1." +
                        " See the documentation of class DGDiffusionIntegrator.")
    parser.add_argument("-k", "--kappa",
                        action='store', default=-1.0, type=float,
                        help="One of the two DG penalty parameters, should be positive." +
                        " Negative values are replaced with (order+1)^2.")
    parser.add_argument("-rs", "--refine-serial",
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly in serial.")
    parser.add_argument("-rp", "--refine-parallel",
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly in parallel.")
    parser.add_argument("-mat", "--material-value",
                        action='store', default=1.0, type=float,
                        help="Constant value for material coefficient " +
                        "in the Laplace operator.")
    parser.add_argument("-dbc", "--dirichlet-value",
                        action='store', default=0.0, type=float,
                        help="Constant value for Dirichlet Boundary Condition.")
    parser.add_argument("-nbc", "--neumann-value",
                        action='store', default=1.0, type=float,
                        help="Constant value for Neumann Boundary Condition.")
    parser.add_argument("-rbc-a", "--robin-a-value",
                        action='store', default=1.0, type=float,
                        help="Constant 'a' value for Robin Boundary Condition. " +
                        "du/dn + a * u = b.")
    parser.add_argument("-rbc-b", "--robin-b-value",
                        action='store', default=1.0, type=float,
                        help="Constant 'b' value for Robin Boundary Condition: " +
                        "du/dn + a * u = b.")
    parser.add_argument("-a", "--radius",
                        action='store', default=0.2, type=float,
                        help="Radius of holes in the mesh.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true', default=True,
                        help='Enable GLVis visualization')
    parser.add_argument('-no-vis', '--no_visualization',
                        action='store_true', default=False,
                        help='Disable GLVis visualization')

    args = parser.parse_args()

    h1 = args.continuous
    if args.discontinuous:
        h1 = False
        args.continuous = False
    vis = True
    if args.no_visualization:
        vis = False
        args.visualization = False

    if myid == 0:
        parser.print_options(args)
    run(order=args.order,
        h1=h1,
        sigma=args.sigma,
        kappa=args.kappa,
        ser_ref_levels=args.refine_serial,
        par_ref_levels=args.refine_parallel,
        mat_val=args.material_value,
        dbc_val=args.dirichlet_value,
        nbc_val=args.neumann_value,
        rbc_a_val=args.robin_a_value,
        rbc_b_val=args.robin_b_value,
        a_=args.radius,
        visualization=vis)
