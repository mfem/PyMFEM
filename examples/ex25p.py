
'''
   MFEM example 25
      See c++ version in the MFEM library for more detail

      This example also demonstrates how to use JITed function in
      more genral form.

      In particular it shows.
         * mfem.jit.vector decorator as a function call
         * defining constant parameters in functions
         * calling JITed function from JITed coefficient
     
'''
from numba import jit, types, carray
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi
import scipy.special
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)

prob = ''


def print0(*args, **kwargs):
    if myid == 0:
        print(*args, **kwargs)
    else:
        pass


def run(meshfile="",
        order=1,
        ref_levels=0,
        par_ref_levels=0,
        visualization=1,
        herm_conv=True,
        device_config='cpu',
        pa=False):

    # 2. Enable hardware devices such as GPUs, and programming models such as
    #    CUDA, OCCA, RAJA and OpenMP based on command line options.
    device = mfem.Device(device_config)
    if myid == 0:
        device.Print()

    # 3. Setup the mesh
    if meshfile == '':
        exact_known = True
        if prob == "beam":
            meshfile = "beam-hex.mesh"
        elif prob == "disc":
            meshfile = "square-disc.mesh"
        elif prob == "lshape":
            meshfile = "l-shape.mesh"
        elif prob == "fichera":
            meshfile = "fichera.mesh"
        else:
            meshfile = "inline-quad.mesh"
            exact_known = False
    else:
        exact_known = True

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', meshfile))

    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    # Setup PML length
    length = np.zeros((dim, 2))

    # 4. Setup the Cartesian PML region.
    if prob == "disc":
        length[:] = 0.2
    elif prob == "lshape":
        length[0, 0] = 0.1
        length[1, 0] = 0.1
    elif prob == "fichera":
        length[0, 1] = 0.5
        length[1, 1] = 0.5
        length[2, 1] = 0.5
    elif prob == "beam":
        length[0, 1] = 2.0
    else:
        length[:] = 0.25

    pml = CartesianPML(mesh, length)
    comp_domain_bdr = pml.comp_dom_bdr
    domain_bdr = pml.dom_bdr

    # 5. Refine the mesh to increase the resolution.
    for l in range(ref_levels):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    for l in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 6. Set element attributes in order to distinguish elements in the
    #    PML region
    pml.SetAttributes(pmesh)

    # 7. Define a finite element space on the mesh. Here we use the Nedelec
    #    finite elements of the specified order.
    fec = mfem.ND_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)

    size = fespace.GlobalTrueVSize()

    print0("Number of finite element unknowns: " + str(size))

    # 8. Determine the list of true essential boundary dofs. In this example,
    #    the boundary conditions are defined based on the specific mesh and the
    #    problem type.

    battrs = pmesh.GetBdrAttributeArray()

    if len(battrs) > 0:
        if prob == "lshape" or prob == "fichera":
            ess_bdr0 = [0]*np.max(battrs)

            for j in range(pmesh.GetNBE()):
                bdrgeom = pmesh.GetBdrElementBaseGeometry(j)
                tr = pmesh.GetBdrElementTransformation(j)
                center = tr.Transform(mfem.Geometries.GetCenter(bdrgeom))

                k = pmesh.GetBdrAttribute(j)
                if prob == "lshape":
                    if (center[0] == 1.0 or center[0] == 0.5 or
                            center[1] == 0.5):
                        ess_bdr0[k - 1] = 1
                else:  # prob == "fichera"
                    if (center[0] == -1.0 or center[0] == 0.0 or
                            center[1] == 0.0 or center[2] == 0.0):
                        ess_bdr0[k - 1] = 1
        else:
            ess_bdr0 = [1]*np.max(battrs)

        ess_bdr = mfem.intArray(ess_bdr0)
    else:
        ess_bdr = mfem.intArray()

    ess_tdof_list = mfem.intArray()
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    # 9. Setup Complex Operator convention
    conv = herm_conv if mfem.ComplexOperator.HERMITIAN else mfem.ComplexOperator.BLOCK_SYMMETRIC

    # 10. Set up the linear form b(.) which corresponds to the right-hand side of
    #     the FEM linear system.

    #     constant parameters used in JITed function can be given by params keyword
    params = {"comp_domain_bdr": comp_domain_bdr,
              "dim": dim,
              "omega": omega,
              "epsilon": epsilon,
              "prob": prob,
              "mu": mu}
    f = mfem.jit.vector(shape=(dim, ), params=params)(source)
    b = mfem.ParComplexLinearForm(fespace, conv)
    if prob == "general":
        b.AddDomainIntegrator(None, mfem.VectorFEDomainLFIntegrator(f))

    b.Assign(0.0)
    b.Assemble()

    # 11. Define the solution vector x as a complex finite element grid function
    #     corresponding to fespace.
    x = mfem.ParComplexGridFunction(fespace)
    x.Assign(0.0)

    #     Here, we JIT compile a fuction using mfem.jit.func and use it later
    #     in a coefficeint.
    sig = types.complex128[:](types.double[:])
    exact_solution = mfem.jit.func(sig, params=params)(maxwell_solution)

    params = {'comp_domain_bdr': comp_domain_bdr,
              'exact_solution': exact_solution}
    E_data = mfem.jit.vector(shape=(dim,),
                             complex=True,
                             params=params)(E_bdr_data)
    x.ProjectBdrCoefficientTangent(E_data.real, E_data.imag, ess_bdr)

    # 12. Set up the sesquilinear form a(.,.)
    #
    #     In Comp
    #     Domain:   1/mu (Curl E, Curl F) - omega^2 * epsilon (E,F)
    #
    #     In PML:   1/mu (1/det(J) J^T J Curl E, Curl F)
    #               - omega^2 * epsilon (det(J) * (J^T J)^-1 * E, F)
    #
    #     where J denotes the Jacobian Matrix of the PML Stretching function

    attrs = pmesh.GetAttributeArray()
    if len(attrs) > 0:
        attr = [0]*np.max(attrs)
        attrPML = [0]*np.max(attrs)

        attr[0] = 1
        if max(attrs) > 1:
            attrPML[1] = 1

    muinv = mfem.ConstantCoefficient(1/mu)
    omeg = mfem.ConstantCoefficient(-omega**2 * epsilon)
    attr = mfem.intArray(attr)
    attrPML = mfem.intArray(attrPML)
    restr_muinv = mfem.RestrictedCoefficient(muinv, attr)
    restr_omeg = mfem.RestrictedCoefficient(omeg, attr)

    # Integrators inside the computational domain (excluding the PML region)
    a = mfem.ParSesquilinearForm(fespace, conv)
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(restr_muinv), None)
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(restr_omeg), None)

    cdim = 1 if dim == 2 else dim

    # JIT compiles all functions first. params defines local variables
    # inside the JITed function.
    params = {"StretchFunction": pml.StretchFunction}
    detJ_inv_JT_J = mfem.jit.vector(shape=(dim,),
                                    params=params,
                                    complex=True)(detJ_inv_JT_J_f)
    detJ_inv_JT_J_abs = mfem.jit.vector(shape=(dim,),
                                        params=params)(detJ_inv_JT_J_abs_f)

    detJ_JT_J_inv = mfem.jit.vector(params=params,
                                    shape=(dim,),
                                    complex=True)(detJ_JT_J_inv_f)
    detJ_JT_J_inv_abs = mfem.jit.vector(shape=(dim,),
                                        params=params)(detJ_JT_J_inv_abs_f)

    def dm1(x, diag):
        return (1/mu)*diag

    c1 = mfem.jit.vector(shape=(cdim, ),
                         dependency=(detJ_inv_JT_J,), complex=True)(dm1)
    c1_Re = c1.real
    c1_Im = c1.imag
    restr_c1_Re = mfem.VectorRestrictedCoefficient(c1_Re, attrPML)
    restr_c1_Im = mfem.VectorRestrictedCoefficient(c1_Im, attrPML)

    def dm2(x, diag):
        return (-omega**2 * epsilon)*diag

    c2 = mfem.jit.vector(shape=(dim,),
                         dependency=(detJ_JT_J_inv,), complex=True)(dm2)
    c2_Re = c2.real
    c2_Im = c2.imag
    restr_c2_Re = mfem.VectorRestrictedCoefficient(c2_Re, attrPML)
    restr_c2_Im = mfem.VectorRestrictedCoefficient(c2_Im, attrPML)

    # Integrators inside the PML region
    a.AddDomainIntegrator(mfem.CurlCurlIntegrator(restr_c1_Re),
                          mfem.CurlCurlIntegrator(restr_c1_Im))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(restr_c2_Re),
                          mfem.VectorFEMassIntegrator(restr_c2_Im))

    # 13. Assemble the bilinear form and the corresponding linear system,
    #     applying any necessary transformations such as: assembly, eliminating
    #     boundary conditions, applying conforming constraints for
    #     non-conforming AMR, etc.
    if pa:
        a.SetAssemblyLevel(mfem.AssemblyLevel.PARTIAL)
    a.Assemble(0)

    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

    # 14a. Set up the Bilinear form a(.,.) for the preconditioner
    #
    #    In Comp
    #    Domain:   1/mu (Curl E, Curl F) + omega^2 * epsilon (E,F)
    #
    #    In PML:   1/mu (abs(1/det(J) J^T J) Curl E, Curl F)
    #              + omega^2 * epsilon (abs(det(J) * (J^T J)^-1) * E, F)

    umf_solver = False
    if pa or not umf_solver:
        absomeg = mfem.ConstantCoefficient(omega**2 * epsilon)
        restr_absomeg = mfem.RestrictedCoefficient(absomeg, attr)

        prec = mfem.ParBilinearForm(fespace)
        prec.AddDomainIntegrator(mfem.CurlCurlIntegrator(restr_muinv))
        prec.AddDomainIntegrator(mfem.VectorFEMassIntegrator(restr_absomeg))

        c1_abs = mfem.jit.vector(shape=(cdim,),
                                 dependency=(detJ_inv_JT_J_abs,))(dm1)
        restr_c1_abs = mfem.VectorRestrictedCoefficient(c1_abs, attrPML)

        def dm3(x, diag):
            return (omega**2 * epsilon)*diag

        c2_abs = mfem.jit.vector(shape=(dim,),
                                 dependency=(detJ_JT_J_inv_abs,))(dm3)
        restr_c2_abs = mfem.VectorRestrictedCoefficient(c2_abs, attrPML)

        prec.AddDomainIntegrator(mfem.CurlCurlIntegrator(restr_c1_abs))
        prec.AddDomainIntegrator(mfem.VectorFEMassIntegrator(restr_c2_abs))

        if pa:
            prec.SetAssemblyLevel(mfem.AssemblyLevel.PARTIAL)
        prec.Assemble()

        # 14b. Define and apply a GMRES solver for AU=B with a block diagonal
        #      preconditioner based on the Gauss-Seidel or Jacobi sparse smoother.
        offsets = intArray([0, fespace.GetTrueVSize(), fespace.GetTrueVSize()])
        offsets.PartialSum()

        s = -1.0 if conv == mfem.ComplexOperator.HERMITIAN else 1.0
        if pa:
            # Jacobi Smoother
            d00 = mfem.OperatorJacobiSmoother(prec, ess_tdof_list)
            d11 = mfem.ScaledOperator(d00, s)
            pc_r = d00
            pc_i = d11
        else:
            PCOpAh = mfem.OperatorPtr()
            prec.SetDiagonalPolicy(mfem.Operator.DIAG_ONE)
            prec.FormSystemMatrix(ess_tdof_list, PCOpAh)

            # Gauss-Seidel Smoother
            ams00 = mfem.HypreAMS(PCOpAh.AsHypreParMatrix(), fespace)
            ams11 = mfem.ScaledOperator(ams00, s)
            pc_r = ams00
            pc_i = ams11

        BlockDP = mfem.BlockDiagonalPreconditioner(offsets)
        BlockDP.SetDiagonalBlock(0, pc_r)
        BlockDP.SetDiagonalBlock(1, pc_i)

        gmres = mfem.GMRESSolver(MPI.COMM_WORLD)
        gmres.SetPrintLevel(1)
        gmres.SetKDim(200)
        gmres.SetMaxIter(5000 if pa else 2000)
        # gmres.SetMaxIter(1)
        gmres.SetRelTol(1e-5)
        gmres.SetAbsTol(0.0)
        gmres.SetOperator(A.Ptr())
        gmres.SetPreconditioner(BlockDP)
        gmres.Mult(B, X)

    # 15. Recover the solution as a finite element grid function and compute the
    #     errors if the exact solution is known.
    a.RecoverFEMSolution(X, b, x)

    # If exact is known compute the error
    if exact_known:
        E_ex = mfem.jit.vector(shape=(dim,),
                               complex=True)(maxwell_solution)
        E_ex_Re = E_ex.real
        E_ex_Im = E_ex.imag

        order_quad = max([2, 2 * order + 1])

        birs = [mfem.IntRules.Get(i, order_quad)
                for i in range(mfem.Geometry.NumGeom)]

        L2Error_Re = x.real().ComputeL2Error(E_ex_Re, birs, pml.elems)
        L2Error_Im = x.imag().ComputeL2Error(E_ex_Im, birs, pml.elems)

        x_gf0 = mfem.ParComplexGridFunction(fespace)
        x_gf0.Assign(0.0)
        norm_E_Re = x_gf0.real().ComputeL2Error(E_ex_Re, birs, pml.elems)
        norm_E_Im = x_gf0.imag().ComputeL2Error(E_ex_Im, birs, pml.elems)

        print0("")
        print0(" Relative Error (Re part): || E_h - E || / ||E|| = " +
               "{:g}".format(L2Error_Re / norm_E_Re))
        print0(" Relative Error (Im part): || E_h - E || / ||E|| = " +
               "{:g}".format(L2Error_Im / norm_E_Im))
        print0(" Total Error : " +
               "{:g}".format(sqrt(L2Error_Re*L2Error_Re + L2Error_Im*L2Error_Im)))
        print0("")

    pmesh.Print("mesh"+smyid, 8)
    x.real().Save("ex25p-sol_r"+smyid, 8)
    x.imag().Save("ex25p-sol_i"+smyid, 8)

    if visualization:
        keys = "keys macF\n" if dim == 3 else "keys amrRljcUUuu\n"
        if prob == "beam" and dim == 3:
            keys = "keys macFFiYYYYYYYYYYYYYYYYYY\n"
        if prob == "beam" and dim == 2:
            keys = "keys amrRljcUUuuu\n"

        sol_sock_re = mfem.socketstream("localhost", 19916)
        sol_sock_re.precision(8)
        sol_sock_re << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_re << "solution\n" << pmesh << x.real() << keys
        sol_sock_re << "window_title 'Soluiton real part'"
        sol_sock_re.flush()

        sol_sock_im = mfem.socketstream("localhost", 19916)
        sol_sock_im.precision(8)
        sol_sock_im << "parallel " << num_procs << " " << myid << "\n"
        sol_sock_im << "solution\n" << pmesh << x.imag() << keys
        sol_sock_im << "window_title 'Soluiton imag part'"
        sol_sock_im.flush()

        x_t = mfem.ParGridFunction(fespace)
        x_t.Assign(x.real())

        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock << "solution\n" << pmesh << x_t << keys << "autoscale off\n"
        sol_sock << "window_title 'Harmonic Solution (t = 0.0T)'"
        sol_sock << "pause\n"
        sol_sock.flush()

        print0(
            "GLVis visualization paused. Press space (in the GLVis window) to resume it.")
        num_frames = 32
        i = 0

        for i in range(num_frames):
            t = (i % num_frames) / num_frames
            oss = "Harmonic Solution (t = " + str(t) + " T)"
            dd = (cos(2.0 * pi * t)*x.real().GetDataArray() +
                  sin(2.0 * pi * t)*x.imag().GetDataArray())
            # x_t.Assign(mfem.Vector(dd))
            x_t.Assign(dd)
            sol_sock << "parallel " << num_procs << " " << myid << "\n"
            sol_sock << "solution\n" << pmesh << x_t
            sol_sock << "window_title '" << oss << "'"
            sol_sock.flush()


class CartesianPML:
    def __init__(self, mesh, length):
        self.length = length
        self.dim = mesh.Dimension()
        self.SetBoundaries(mesh)

    def SetBoundaries(self, mesh):
        self.comp_dom_bdr = np.zeros((self.dim, 2))
        self.dom_bdr = np.zeros((self.dim, 2))
        pmin, pmax = mesh.GetBoundingBox()
        for i in range(self.dim):
            self.dom_bdr[i, 0] = pmin[i]
            self.dom_bdr[i, 1] = pmax[i]
            self.comp_dom_bdr[i, 0] = self.dom_bdr[i, 0] + self.length[i, 0]
            self.comp_dom_bdr[i, 1] = self.dom_bdr[i, 1] - self.length[i, 1]

    def SetAttributes(self, mesh):  # this mesh is not the same mesh pass to __init__
        # Initialize bdr attributes
        self.elems = mfem.intArray(mesh.GetNE())

        for i in range(mesh.GetNBE()):
            mesh.GetBdrElement(i).SetAttribute(i+1)

        #  Loop through the elements and identify which of them are in the PML
        for i in range(mesh.GetNE()):
            self.elems[i] = 1

            in_pml = False
            el = mesh.GetElement(i)

            #  Initialize attribute
            el.SetAttribute(1)
            vertices = el.GetVerticesArray()
            nrvert = len(vertices)

            # Check if any vertex is in the PML
            for iv in range(nrvert):
                vert_idx = vertices[iv]
                coords = mesh.GetVertexArray(vert_idx)

                for comp in range(self.dim):
                    if (coords[comp] > self.comp_dom_bdr[comp, 1] or
                            coords[comp] < self.comp_dom_bdr[comp, 0]):
                        in_pml = True
                        break

            if in_pml:
                self.elems[i] = 0
                el.SetAttribute(2)

        # construct attribute array in Mesh object
        mesh.SetAttributes()
        self.StretchFunction = self._GenerateStretchFunction()

    def _GenerateStretchFunction(self):
        sig = types.void(types.double[:], types.complex128[:])
        params = {"comp_domain_bdr": self.comp_dom_bdr,
                  "dim": self.dim,
                  "length": self.length,
                  "omega": omega,
                  "epsilon": epsilon,
                  "mu": mu}

        def _StretchFunction(x, dxs):
            zi = 1j

            n = 2.0
            c = 5.0
            k = omega * sqrt(epsilon * mu)

            # Stretch in each direction independently
            for i in range(dim):
                dxs[i] = 1.0
                if x[i] >= comp_domain_bdr[i, 1]:
                    coeff = n * c / k / length[i, 1]**n
                    dxs[i] = (1.0 + zi * coeff *
                              abs((x[i] - comp_domain_bdr[i, 1])**(n-1.0)))
                if x[i] <= comp_domain_bdr[i, 0]:
                    coeff = n * c / k / length[i, 0]**n
                    dxs[i] = (1.0 + zi * coeff *
                              abs((x[i] - comp_domain_bdr[i, 0])**(n-1.0)))

        func = mfem.jit.func(sig, params=params)(_StretchFunction)
        return func

#
#  functions (these are JITed using Numba insdie run())
#


def source(x):
    dim = shape[0]
    out = np.zeros(dim)
    center = np.zeros(dim)
    r = 0
    for i in range(dim):
        center[i] = 0.5 * (comp_domain_bdr[i, 0] + comp_domain_bdr[i, 1])
        r += (x[i] - center[i])**2
        out[i] = 0

    n = 5.0 * omega * sqrt(epsilon * mu) / pi
    coeff = n**2 / pi
    alpha = -n**2 * r

    out[0] = coeff * exp(alpha)
    return out


def maxwell_solution(x):

    jn = scipy.special.jv
    yn = scipy.special.yn

    # Initialize
    E = np.zeros(dim, dtype=np.complex128)

    zi = 1j
    k = omega * sqrt(epsilon * mu)

    if prob == "disc" or prob == "lshape" or prob == "fichera":
        shift = np.zeros(dim)
        if prob == "fichera":
            shift += 1.0
        elif prob == "disc":
            shift -= 0.5
        else:
            shift -= 1.0

        if dim == 2:
            x0 = x[0] + shift[0]
            x1 = x[1] + shift[1]
            r = sqrt(x0 * x0 + x1 * x1)
            beta = k * r

            # Bessel functions
            Ho = jn(0.0, beta) + zi * yn(0, beta)
            Ho_r = -k * (jn(1.0, beta) + zi * yn(1, beta))
            Ho_rr = -k * k * (1.0 / beta *
                              (jn(1., beta) + zi * yn(1, beta)) -
                              (jn(2., beta) + zi * yn(2, beta)))

            # First derivatives
            r_x = x0 / r
            r_y = x1 / r
            r_xy = -(r_x / r) * r_y
            r_xx = (1.0 / r) * (1.0 - r_x * r_x)

            val = 0.25 * zi * Ho
            val_xx = 0.25 * zi * (r_xx * Ho_r + r_x * r_x * Ho_rr)
            val_xy = 0.25 * zi * (r_xy * Ho_r + r_x * r_y * Ho_rr)
            E[0] = zi / k * (k * k * val + val_xx)
            E[1] = zi / k * val_xy
        elif dim == 3:
            x0 = x[0] + shift[0]
            x1 = x[1] + shift[1]
            x2 = x[2] + shift[2]
            r = sqrt(x0 * x0 + x1 * x1 + x2 * x2)

            r_x = x0 / r
            r_y = x1 / r
            r_z = x2 / r

            r_xx = (1.0 / r) * (1.0 - r_x * r_x)
            r_yx = -(r_y / r) * r_x
            r_zx = -(r_z / r) * r_x

            val = exp(zi * k * r) / r
            val_r = val / r * (zi * k * r - 1.0)
            val_rr = val / (r * r) * (-k * k * r * r
                                      - 2.0 * zi * k * r + 2.0)

            val_xx = val_rr * r_x * r_x + val_r * r_xx
            val_yx = val_rr * r_x * r_y + val_r * r_yx
            val_zx = val_rr * r_x * r_z + val_r * r_zx

            alpha = zi * k / 4.0 / pi / k / k

            E[0] = alpha * (k * k * val + val_xx)
            E[1] = alpha * val_yx
            E[2] = alpha * val_zx
        else:
            pass
    elif prob == 'beam':
        # T_10 mode
        if dim == 3:
            k10 = sqrt(k * k - pi*pi)
            E[1] = -zi * k / pi * sin(pi*x[2])*exp(zi * k10 * x[0])
        elif dim == 2:
            E[1] = -zi * k / pi * exp(zi * k * x[0])
        else:
            pass
    else:
        pass
    return E


def E_bdr_data(x):
    dim = shape[0]
    E = np.zeros(dim, dtype=np.complex128)

    in_pml = False

    for i in range(dim):
        # check if in PML
        if ((x[i] - comp_domain_bdr[i, 0]) < 0.0 or
                (x[i] - comp_domain_bdr[i, 1]) > 0.0):
            in_pml = True
            break
    if not in_pml:
        E = exact_solution(x)

    return E


def detJ_JT_J_inv_f(x):
    dim = shape[0]
    D = np.zeros(dim, dtype=np.complex128)
    # print(D) # this print cause memory error....(why?)
    dxs = np.empty(dim, dtype=np.complex128)
    det = complex(1.0)

    StretchFunction(x, dxs)
    for i in range(dim):
        det *= dxs[i]
    for i in range(dim):
        D[i] = (det / (dxs[i]**2))

    return D


def detJ_JT_J_inv_abs_f(x):
    dim = shape[0]
    D = np.zeros(dim, dtype=np.float64)

    dxs = np.empty(dim, dtype=np.complex128)
    det = complex(1.0)
    StretchFunction(x, dxs)
    for i in range(dim):
        det *= dxs[i]
    for i in range(dim):
        D[i] = abs(det / (dxs[i]**2))
    return D


def detJ_inv_JT_J_f(x):
    dim = shape[0]
    D = np.zeros(dim, dtype=np.complex128)

    dxs = np.empty(dim, dtype=np.complex128)
    det = 1.0
    StretchFunction(x, dxs)
    for i in range(dim):
        det *= dxs[i]
    # in the 2D case the coefficient is scalar 1/det(J)
    if dim == 2:
        D[0] = (1.0 / det)
    else:
        for i in range(dim):
            D[i] = (dxs[i]**2 / det)
    return D


def detJ_inv_JT_J_abs_f(x):
    dim = shape[0]
    D = np.zeros(dim, dtype=np.float64)

    dxs = np.empty(dim, dtype=np.complex128)
    det = 1.0
    StretchFunction(x, dxs)
    for i in range(dim):
        det *= dxs[i]
    # in the 2D case the coefficient is scalar 1/det(J)
    if dim == 2:
        D[0] = abs(1.0 / det)
    else:
        for i in range(dim):
            D[i] = abs(dxs[i]**2 / det)
    return D


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex25 (PML)')
    parser.add_argument('-m', '--mesh',
                        default="",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-prob",
                        "--problem-type",
                        action='store', type=int, default=4,
                        help=" 0: beam, 1: disc, 2: lshape, 3: fichera, 4: General")
    parser.add_argument("-rs", "--refinements-serial",
                        action='store', type=int, default=1,
                        help="Number of serial refinements")
    parser.add_argument("-rp", "--refinements-parallel",
                        action='store', type=int, default=2,
                        help="Number of parallel refinements")
    parser.add_argument("-mu", "--permeability",
                        action='store', type=float, default=1.0,
                        help="Permeability of free space (or 1/(spring constant)).")
    parser.add_argument("-eps", "--permittivity",
                        action='store', type=float, default=1.0,
                        help="Permittivity of free space (or mass constant).")
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=5.0,
                        help="Set the frequency for the exact")
    parser.add_argument("-no-herm", "--no-hermitian",
                        action='store_false',
                        default=True,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=True,
                        help='Enable GLVis visualization')
    parser.add_argument("-pa", "--partial-assembly",
                        action='store_true',
                        help="Enable Partial Assembly.")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    probs = {0: "beam", 1: "disc", 2: "lshape", 3: "fichera", 4: "general"}
    globals()["prob"] = probs[args.problem_type]
    globals()["omega"] = 2*pi*args.frequency
    globals()["epsilon"] = args.permittivity
    globals()["mu"] = args.permeability
    run(meshfile=args.mesh,
        order=args.order,
        ref_levels=args.refinements_serial,
        par_ref_levels=args.refinements_parallel,
        visualization=args.visualization,
        herm_conv=args.no_hermitian,
        device_config=args.device,
        pa=args.partial_assembly)
