##
# MFEM Ultraweak DPG Maxwell parallel example
##
# Compile with: make pmaxwell
##
# sample run

# mpirun -np 1 python pmaxwell.py  -o 2 -pref 2 -prob 0
# mpirun -np 4 python pmaxwell.py  -o 3 -sref 0 -pref 3 -rnum 4.8 -sc -prob 0
# mpirun -np 1 python pmaxwell.py  -o 2 -pref 2 -prob 0 -m inline-hex.mesh

# mpirun -np 4 python pmaxwell.py -m ../../data/star.mesh -o 2 -sref 0 -pref 3 -rnum 0.5 -prob 0
# mpirun -np 4 python pmaxwell.py -m ../../data/inline-quad.mesh -o 3 -sref 0 -pref 3 -rnum 4.8 -sc -prob 0
# mpirun -np 4 python pmaxwell.py -m ../../data/inline-hex.mesh -o 2 -sref 0 -pref 1 -rnum 0.8 -sc -prob 0
# mpirun -np 4 python pmaxwell.py -m ../../data/inline-quad.mesh -o 3 -sref 1 -pref 3 -rnum 4.8 -sc -prob 2
# mpirun -np 4 python pmaxwell.py -o 3 -sref 1 -pref 2 -rnum 11.8 -sc -prob 3
# mpirun -np 4 python pmaxwell.py -o 3 -sref 1 -pref 2 -rnum 9.8 -sc -prob 4

# AMR run. Note that this is a computationally intensive sample run.
# We recommend trying it on a large machine with more mpi ranks
# mpirun -np 4 pmaxwell -o 3 -sref 0 -pref 15 -prob 1 -theta 0.7 -sc

# Description:
# This example code demonstrates the use of MFEM to define and solve
# the "ultraweak" (UW) DPG formulation for the Maxwell problem

# ∇×(1/μ ∇×E) - ω² ϵ E = Ĵ ,   in Ω
# E×n = E₀ , on ∂Ω

# It solves the following kinds of problems
# 1) Known exact solutions with error convergence rates
# a) A manufactured solution problem where E is a plane beam
# 2) Fichera "microwave" problem
# 3) PML problems
# a) Generic PML problem with point source given by the load
# b) Plane wave scattering from a square
# c) PML problem with a point source prescribed on the boundary

# The DPG UW deals with the First Order System
# i ω μ H + ∇ × E = 0,   in Ω
# -i ω ϵ E + ∇ × H = J,   in Ω
# E × n = E_0, on ∂Ω
# Note: Ĵ = -iωJ

# The ultraweak-DPG formulation is obtained by integration by parts of both
# equations and the introduction of trace unknowns on the mesh skeleton

# in 2D
# E is vector valued and H is scalar.
# (∇ × E, F) = (E, ∇ × F) + < n × E , F>
# or (∇ ⋅ AE , F) = (AE, ∇ F) + < AE ⋅ n, F>
# where A = [0 1; -1 0];

# E ∈ (L²(Ω))² , H ∈ L²(Ω)
# Ê ∈ H^-1/2(Γₕ), Ĥ ∈ H^1/2(Γₕ)
# i ω μ (H,F) + (E, ∇ × F) + < AÊ, F > = 0,      ∀ F ∈ H¹
# -i ω ϵ (E,G) + (H,∇ × G)  + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
# Ê = E₀      on ∂Ω
# -------------------------------------------------------------------------
# |   |       E      |      H      |      Ê       |       Ĥ      |  RHS    |
# -------------------------------------------------------------------------
# | F |  (E,∇ × F)   | i ω μ (H,F) |   < Ê, F >   |              |         |
# |   |              |             |              |              |         |
# | G | -i ω ϵ (E,G) |  (H,∇ × G)  |              | < Ĥ, G × n > |  (J,G)  |
# where (F,G) ∈  H¹ × H(curl,Ω)

# in 3D
# E,H ∈ (L^2(Ω))³
# Ê ∈ H_0^1/2(Ω)(curl, Γₕ), Ĥ ∈ H^-1/2(curl, Γₕ)
# i ω μ (H,F) + (E,∇ × F) + < Ê, F × n > = 0,      ∀ F ∈ H(curl,Ω)
# -i ω ϵ (E,G) + (H,∇ × G) + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
# Ê × n = E₀      on ∂Ω
# -------------------------------------------------------------------------
# |   |       E      |      H      |      Ê       |       Ĥ      |  RHS    |
# -------------------------------------------------------------------------
# | F |  (E,∇ × F)   | i ω μ (H,F) | < n × Ê, F > |              |         |
# |   |              |             |              |              |         |
# | G | -i ω ϵ (E,G) |  (H,∇ × G)  |              | < n × Ĥ, G > |  (J,G)  |
# where (F,G) ∈  H(curl,Ω) × H(curl,Ω)

# Here we use the "Adjoint Graph" norm on the test space i.e.,
# ||(F,G)||²ᵥ  = ||A^*(F,G)||² + ||(F,G)||² where A is the
# maxwell operator defined by (1)

# The PML formulation is

# ∇×(1/μ α ∇×E) - ω² ϵ β E = Ĵ ,   in Ω
# E×n = E₀ , on ∂Ω

# where α = |J|⁻¹ Jᵀ J (in 2D it's the scalar |J|⁻¹),
# β = |J| J⁻¹ J⁻ᵀ, J is the Jacobian of the stretching map
# and |J| its determinant.

# The first order system reads
# i ω μ α⁻¹ H + ∇ × E = 0,   in Ω
# -i ω ϵ β E + ∇ × H = J,   in Ω
# E × n = E₀,  on ∂Ω

# and the ultraweak formulation is

# in 2D
# E ∈ (L²(Ω))² , H ∈ L²(Ω)
# Ê ∈ H^-1/2(Ω)(Γₕ), Ĥ ∈ H^1/2(Γₕ)
# i ω μ (α⁻¹ H,F) + (E, ∇ × F) + < AÊ, F > = 0,          ∀ F ∈ H¹
# -i ω ϵ (β E,G)   + (H,∇ × G)  + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
# Ê = E₀     on ∂Ω
# ---------------------------------------------------------------------------------
# |   |       E        |        H         |      Ê       |       Ĥ      |  RHS    |
# ---------------------------------------------------------------------------------
# | F |  (E,∇ × F)     | i ω μ (α⁻¹ H,F)  |   < Ê, F >   |              |         |
# |   |                |                  |              |              |         |
# | G | -i ω ϵ (β E,G) |    (H,∇ × G)     |              | < Ĥ, G × n > |  (J,G)  |

# where (F,G) ∈  H¹ × H(curl,Ω)

##
# in 3D
# E,H ∈ (L^2(Ω))³
# Ê ∈ H_0^1/2(Ω)(curl, Γ_h), Ĥ ∈ H^-1/2(curl, Γₕ)
# i ω μ (α⁻¹ H,F) + (E,∇ × F) + < Ê, F × n > = 0,      ∀ F ∈ H(curl,Ω)
# -i ω ϵ (β E,G)    + (H,∇ × G) + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
# Ê × n = E_0     on ∂Ω
# -------------------------------------------------------------------------------
# |   |       E      |      H           |      Ê       |       Ĥ      |  RHS    |
# -------------------------------------------------------------------------------
# | F |  ( E,∇ × F)  | i ω μ (α⁻¹ H,F)  | < n × Ê, F > |              |         |
# |   |              |                  |              |              |         |
# | G | -iωϵ (β E,G) |   (H,∇ × G)      |              | < n × Ĥ, G > |  (J,G)  |
# where (F,G) ∈  H(curl,Ω) × H(curl,Ω)

# For more information see https:##doi.org/10.1016/j.camwa.2021.01.017
from mpi4py import MPI
from numba import njit, void, int32, int64, float64, complex128, types
from mfem.common.bessel import yv as yn
from mfem.common.bessel import jv as jn
import os
from os.path import expanduser, join, dirname

import numpy as np
from numpy import pi, exp

import mfem.par as mfem
mdpg = mfem.dpg


num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)

visport = 19916


def VisualizeField(sock, vishost, visport,
                   gf, title,
                   x=0, y=0, w=400, h=400, keys='', vec=False):

    pmesh = gf.ParFESpace().GetParMesh()

    newly_opened = False

    while True:
        if sock is None:
            sock = mfem.socketstream("localhost", visport)
            sock.precision(8)
            newly_opened = True

        sock << "parallel " << str(num_procs) << " " << str(myid) << "\n"
        sock << "solution\n"
        sock << pmesh << gf

        if newly_opened:
            sock << "window_title '" << title << "'\n"
            sock << "window_geometry "
            sock << x << " " << y << " " << w << " " << h << "\n"

            if keys != '':
                sock << "keys " << keys << "\n"
            else:
                sock << "keys maaAc\n"
            if vec:
                sock << "vvv"
            sock.flush()
        if sock.is_open():
            break

    return sock


prob_type = ("plane_wave",
             "fichera_oven",
             "pml_general",
             "pml_plane_wave_scatter",
             "pml_pointsource")


class Functions:
    def __init__(self, dim, dimc, omega, epsilon, mu, iprob):

        prob = prob_type[iprob]

        @njit(complex128[:](float64[:]))
        def maxwell_solution(x):
            E = np.zeros(dim, dtype=np.complex128)

            if prob == "plane_wave":
                E[0] = exp(1j * omega * np.sum(x))

            elif prob == "pml_plane_wave_scatter":
                E[1] = exp(1j * omega * x[0])

            elif prob == "fichera_oven":
                if abs(x[2] - 3.0) < 1e-10:
                    E[0] = np.sin(pi*x[1])
            elif prob == "pml_pointsource":
                k = omega * np.sqrt(epsilon * mu)
                shift = np.array([-0.5]*dim)

                if dim == 2:
                    x0 = x[0] + shift[0]
                    x1 = x[1] + shift[1]
                    r = np.sqrt(x0 * x0 + x1 * x1)
                    beta = k * r

                    # Bessel functions
                    Ho = jn(0, beta).real + 1j * yn(0, beta).real
                    Ho_r = -k * jn(1, beta).real + 1j * yn(1, beta).real
                    Ho_rr = -k * k * (1. / beta *
                                      (jn(1, beta).real + 1j * yn(1, beta).real) -
                                      (jn(2, beta).real + 1j * yn(2, beta).real))

                    # First derivatives
                    r_x = x0 / r
                    r_y = x1 / r
                    r_xy = -(r_x / r) * r_y
                    r_xx = (1.0 / r) * (1.0 - r_x * r_x)

                    val = 0.25 * 1j * Ho
                    val_xx = 0.25 * 1j * (r_xx * Ho_r + r_x * r_x * Ho_rr)
                    val_xy = 0.25 * 1j * (r_xy * Ho_r + r_x * r_y * Ho_rr)
                    E[0] = 1j / k * (k * k * val + val_xx)
                    E[1] = 1j / k * val_xy
                else:
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

                    val = exp(1j * k * r) / r
                    val_r = val / r * (1j * k * r - 1.)
                    val_rr = val / (r * r) * (-k * k * r * r
                                              - 2. * 1j * k * r + 2.)

                    val_xx = val_rr * r_x * r_x + val_r * r_xx
                    val_yx = val_rr * r_x * r_y + val_r * r_yx
                    val_zx = val_rr * r_x * r_z + val_r * r_zx
                    alpha = 1j * k / 4. / pi / k / k
                    E[0] = alpha * (k * k * val + val_xx)
                    E[1] = alpha * val_yx
                    E[2] = alpha * val_zx
            else:
                assert False, "should not come here"
            return E

        @njit(complex128[:](float64[:]))
        def maxwell_solution_curl(x):
            curlE = np.zeros(dimc, dtype=np.complex128)
            if prob == "plane_wave":
                pw = exp(1j * omega * np.sum(x))
                if dim == 3:
                    curlE[0] = 0.0
                    curlE[1] = 1j * omega * pw
                    curlE[2] = -1j * omega * pw
                else:
                    curlE[0] = -1j * omega * pw
            elif prob == "pml_plane_wave_scatter":
                pw = exp(1j * omega * (x[0]))
                curlE[0] = 1j * omega * pw
            else:
                assert False, "should not come here"
            return curlE

        @njit(complex128[:](float64[:]))
        def maxwell_solution_curlcurl(x):
            curlcurlE = np.zeros(dim, dtype=np.complex128)
            if prob == "plane_wave":
                pw = exp(1j * omega * np.sum(x))
                if dim == 3:
                    curlcurlE[0] = 2. * omega * omega * pw
                    curlcurlE[1] = - omega * omega * pw
                    curlcurlE[2] = - omega * omega * pw
                else:
                    curlcurlE[0] = omega * omega * pw
                    curlcurlE[1] = -omega * omega * pw
            elif prob == "pml_plane_wave_scatter":
                pw = np.exp(1j * omega * x[0])
                curlcurlE[1] = omega * omega * pw
            else:
                assert False, "should not come here"
            return curlcurlE

        self.maxwell_solution = maxwell_solution
        self.maxwell_solution_curl = maxwell_solution_curl
        self.maxwell_solution_curlcurl = maxwell_solution_curlcurl

        @njit(complex128[:](float64[:]))
        def curlH_exact(x):
            # ∇ × H = ∇ × ∇ × E / ω μ
            curlcurlE = maxwell_solution_curlcurl(x)
            return 1j*curlcurlE / (omega * mu)

        @njit(complex128[:](float64[:]))
        def H_exact(x):
            H = maxwell_solution_curl(x)*1j/omega/mu
            return H

        @njit(complex128[:](float64[:]))
        def hatE_exact(x):
            if dim == 3:
                E = maxwell_solution(x)
                return E
            else:
                ret = np.zeros(2, dtype=np.complex128)
                E = maxwell_solution(x)
                ret[0] = E[1]
                ret[1] = -E[0]
            return ret

        @mfem.jit.vector(vdim=dim)
        def E_exact_r(x):
            E = maxwell_solution(x)
            return E.real

        @mfem.jit.vector(vdim=dim)
        def E_exact_i(x):
            E = maxwell_solution(x)
            return E.imag

        @mfem.jit.vector(vdim=dim)
        def H_exact_r(x):
            return H_exact(x).real

        @mfem.jit.vector(vdim=dim)
        def H_exact_i(x):
            return H_exact(x).imag

        @mfem.jit.vector(vdim=dim)
        def hatE_exact_r(x):
            ret = hatE_exact(x)
            return ret.real

        @mfem.jit.vector(vdim=dim)
        def hatE_exact_i(x):
            ret = hatE_exact(x)
            return ret.imag

        self.E_exact_r = E_exact_r
        self.E_exact_i = E_exact_i
        self.H_exact_r = H_exact_r
        self.H_exact_i = H_exact_i
        self.hatE_exact_r = hatE_exact_r
        self.hatE_exact_i = hatE_exact_i

        # J = -i ω ϵ E + ∇ × H
        # J_r + iJ_i = -i ω ϵ (E_r + i E_i) + ∇ × (H_r + i H_i)
        @mfem.jit.vector(vdim=dim)
        def rhs_func_r(x):
            # J_r = ω ϵ E_i + ∇ × H_r
            E = maxwell_solution(x)
            curlH = curlH_exact(x)
            return omega * epsilon * E.imag + curlH.real

        @mfem.jit.vector(vdim=dim)
        def rhs_func_i(x):
            # J_r = ω ϵ E_i + ∇ × H_r
            E = maxwell_solution(x)
            curlH = curlH_exact(x)
            return -omega * epsilon * E.real + curlH.imag

        self.rhs_func_r = rhs_func_r
        self.rhs_func_i = rhs_func_i

        @mfem.jit.vector(vdim=dim)
        def source_function(x):
            center = np.zeros(dim)+0.5
            r = 0.0
            for i in range(dim):
                r += (x[i] - center[i])**2

            n = 5.0 * omega * np.sqrt(epsilon * mu) / pi
            coeff = pow(n, 2) / pi
            alpha = -pow(n, 2) * r
            f = np.zeros(dim)
            f[0] = -omega * coeff * exp(alpha)/omega
            return f

        self.source_function = source_function


def run(meshfile='',
        order=1,
        delta_order=1,
        prob=0,
        sr=0,
        pr=1,
        epsilon=1.0,
        mu=1.0,
        rnum=1.0,
        theta=0.0,
        static_cond=False,
        visualization=False):

    omega = 2.*pi*rnum
    with_pml = False
    exact_known = False

    if prob == 0:
        exact_known = True
        mesh_file = expanduser(
            join(dirname(__file__), '..', '..', 'data', meshfile))

    elif prob == 1:
        meshfile = "meshes/fichera-waveguide.mesh"
        omega = 5.0
        rnum = omega/(2.*pi)
        mesh_file = expanduser(
            join(dirname(__file__),  meshfile))

    elif prob == 2:
        with_pml = True
        mesh_file = expanduser(
            join(dirname(__file__), '..', '..', 'data', meshfile))

    else:
        with_pml = True
        meshfile = "meshes/scatter.mesh"
        mesh_file = expanduser(
            join(dirname(__file__),  meshfile))

    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()
    assert dim > 1, "Dimension = 1 is not supported in this example"

    dimc = 3 if dim == 3 else 1
    for i in range(sr):
        mesh.UniformRefinement()
    mesh.EnsureNCMesh(False)

    pml = None
    if with_pml:
        length = mfem.doubleArray2D(dim, 2)
        length.Assign(0.25)
        pml = mdpg.CartesianPML(mesh, length)
        pml.SetOmega(omega)
        pml.SetEpsilonAndMu(epsilon, mu)

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh

    attr = mfem.intArray()
    attrPML = mfem.intArray()
    if with_pml:
        pml.SetAttributes(pmesh, attr, attrPML)

    # Define spaces
    TrialSpace = {"E_space": 0,
                  "H_space": 1,
                  "hatE_space": 2,
                  "hatH_space": 3, }
    TestSpace = {"F_space": 0,
                 "G_space": 1, }

    # Vector L2 L2 space for E
    E_fec = mfem.L2_FECollection(order-1, dim)
    E_fes = mfem.ParFiniteElementSpace(pmesh, E_fec, dim)

    # Vector L2 L2 space for H
    H_fec = mfem.L2_FECollection(order-1, dim)
    H_fes = mfem.ParFiniteElementSpace(pmesh, H_fec, dimc)

    # H^-1/2 (curl) space for Ê
    test_order = order + delta_order

    if dim == 3:
        hatE_fec = mfem.ND_Trace_FECollection(order, dim)
        hatH_fec = mfem.ND_Trace_FECollection(order, dim)
        F_fec = mfem.ND_FECollection(test_order, dim)
    else:
        hatE_fec = mfem.RT_Trace_FECollection(order-1, dim)
        hatH_fec = mfem.H1_Trace_FECollection(order, dim)
        F_fec = mfem.H1_FECollection(test_order, dim)

    hatE_fes = mfem.ParFiniteElementSpace(pmesh, hatE_fec)
    hatH_fes = mfem.ParFiniteElementSpace(pmesh, hatH_fec)
    G_fec = mfem.ND_FECollection(test_order, dim)

    trial_fes = mfem.ParFiniteElementSpaceArray()
    test_fec = mfem.FiniteElementCollectionArray()
    trial_fes.Append(E_fes)
    trial_fes.Append(H_fes)
    trial_fes.Append(hatE_fes)
    trial_fes.Append(hatH_fes)
    test_fec.Append(F_fec)
    test_fec.Append(G_fec)

    # Bilinear form coefficients
    one = mfem.ConstantCoefficient(1.0)
    eps2omeg2 = mfem.ConstantCoefficient(epsilon*epsilon*omega*omega)
    mu2omeg2 = mfem.ConstantCoefficient(mu*mu*omega*omega)
    muomeg = mfem.ConstantCoefficient(mu*omega)
    negepsomeg = mfem.ConstantCoefficient(-epsilon*omega)
    epsomeg = mfem.ConstantCoefficient(epsilon*omega)
    negmuomeg = mfem.ConstantCoefficient(-mu*omega)

    # for the 2D case
    rot_mat = mfem.DenseMatrix(np.array([[0, 1.], [-1, 0]]))
    rot = mfem.MatrixConstantCoefficient(rot_mat)
    epsrot = mfem.ScalarMatrixProductCoefficient(epsomeg, rot)
    negepsrot = mfem.ScalarMatrixProductCoefficient(negepsomeg, rot)

    if with_pml:
        epsomeg_cf = mfem.RestrictedCoefficient(epsomeg, attr)
        negepsomeg_cf = mfem.RestrictedCoefficient(negepsomeg, attr)
        eps2omeg2_cf = mfem.RestrictedCoefficient(eps2omeg2, attr)
        muomeg_cf = mfem.RestrictedCoefficient(muomeg, attr)
        negmuomeg_cf = mfem.RestrictedCoefficient(negmuomeg, attr)
        mu2omeg2_cf = mfem.RestrictedCoefficient(mu2omeg2, attr)
        epsrot_cf = mfem.MatrixRestrictedCoefficient(epsrot, attr)
        negepsrot_cf = mfem.MatrixRestrictedCoefficient(negepsrot, attr)
    else:
        epsomeg_cf = epsomeg
        negepsomeg_cf = negepsomeg
        eps2omeg2_cf = eps2omeg2
        muomeg_cf = muomeg
        negmuomeg_cf = negmuomeg
        mu2omeg2_cf = mu2omeg2
        epsrot_cf = epsrot
        negepsrot_cf = negepsrot

    detJ_r = mdpg.PmlCoefficient(mdpg.detJ_r_function, pml)
    detJ_i = mdpg.PmlCoefficient(mdpg.detJ_i_function, pml)
    abs_detJ_2 = mdpg.PmlCoefficient(mdpg.abs_detJ_2_function, pml)
    detJ_Jt_J_inv_r = mdpg.PmlMatrixCoefficient(
        dim, mdpg.detJ_Jt_J_inv_r_function, pml)
    detJ_Jt_J_inv_i = mdpg.PmlMatrixCoefficient(
        dim, mdpg.detJ_Jt_J_inv_i_function, pml)
    abs_detJ_Jt_J_inv_2 = mdpg.PmlMatrixCoefficient(
        dim, mdpg.abs_detJ_Jt_J_inv_2_function, pml)
    negmuomeg_detJ_r = mfem.ProductCoefficient(negmuomeg, detJ_r)
    negmuomeg_detJ_i = mfem.ProductCoefficient(negmuomeg, detJ_i)
    muomeg_detJ_r = mfem.ProductCoefficient(muomeg, detJ_r)
    mu2omeg2_detJ_2 = mfem.ProductCoefficient(mu2omeg2, abs_detJ_2)
    epsomeg_detJ_Jt_J_inv_i = mfem.ScalarMatrixProductCoefficient(
        epsomeg, detJ_Jt_J_inv_i)
    epsomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(
        epsomeg, detJ_Jt_J_inv_r)
    negepsomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(
        negepsomeg, detJ_Jt_J_inv_r)
    muomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(
        muomeg, detJ_Jt_J_inv_r)
    negmuomeg_detJ_Jt_J_inv_i = mfem.ScalarMatrixProductCoefficient(
        negmuomeg, detJ_Jt_J_inv_i)
    negmuomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(
        negmuomeg, detJ_Jt_J_inv_r)
    mu2omeg2_detJ_Jt_J_inv_2 = mfem.ScalarMatrixProductCoefficient(
        mu2omeg2, abs_detJ_Jt_J_inv_2)
    eps2omeg2_detJ_Jt_J_inv_2 = mfem.ScalarMatrixProductCoefficient(
        eps2omeg2, abs_detJ_Jt_J_inv_2)
    negmuomeg_detJ_r_restr = mfem.RestrictedCoefficient(
        negmuomeg_detJ_r, attrPML)
    negmuomeg_detJ_i_restr = mfem.RestrictedCoefficient(
        negmuomeg_detJ_i, attrPML)
    muomeg_detJ_r_restr = mfem.RestrictedCoefficient(muomeg_detJ_r, attrPML)
    mu2omeg2_detJ_2_restr = mfem.RestrictedCoefficient(
        mu2omeg2_detJ_2, attrPML)
    epsomeg_detJ_Jt_J_inv_i_restr = mfem.MatrixRestrictedCoefficient(
        epsomeg_detJ_Jt_J_inv_i, attrPML)
    epsomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(
        epsomeg_detJ_Jt_J_inv_r, attrPML)
    negepsomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(
        negepsomeg_detJ_Jt_J_inv_r, attrPML)
    muomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(
        muomeg_detJ_Jt_J_inv_r, attrPML)
    negmuomeg_detJ_Jt_J_inv_i_restr = mfem.MatrixRestrictedCoefficient(
        negmuomeg_detJ_Jt_J_inv_i, attrPML)
    negmuomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(
        negmuomeg_detJ_Jt_J_inv_r, attrPML)
    mu2omeg2_detJ_Jt_J_inv_2_restr = mfem.MatrixRestrictedCoefficient(
        mu2omeg2_detJ_Jt_J_inv_2, attrPML)
    eps2omeg2_detJ_Jt_J_inv_2_restr = mfem.MatrixRestrictedCoefficient(
        eps2omeg2_detJ_Jt_J_inv_2, attrPML)

    if with_pml and dim == 2:
        epsomeg_detJ_Jt_J_inv_i_rot = mfem.MatrixProductCoefficient(
            epsomeg_detJ_Jt_J_inv_i, rot)
        epsomeg_detJ_Jt_J_inv_r_rot = mfem.MatrixProductCoefficient(
            epsomeg_detJ_Jt_J_inv_r, rot)
        negepsomeg_detJ_Jt_J_inv_r_rot = mfem.MatrixProductCoefficient(
            negepsomeg_detJ_Jt_J_inv_r, rot)
        epsomeg_detJ_Jt_J_inv_i_rot_restr = mfem.MatrixRestrictedCoefficient(epsomeg_detJ_Jt_J_inv_i_rot,
                                                                             attrPML)
        epsomeg_detJ_Jt_J_inv_r_rot_restr = mfem.MatrixRestrictedCoefficient(epsomeg_detJ_Jt_J_inv_r_rot,
                                                                             attrPML)
        negepsomeg_detJ_Jt_J_inv_r_rot_restr = mfem.MatrixRestrictedCoefficient(negepsomeg_detJ_Jt_J_inv_r_rot,
                                                                                attrPML)

    a = mdpg.ParComplexDPGWeakForm(trial_fes, test_fec)
    a.StoreMatrices()  # needed for AMR

    # (E,∇ × F)
    a.AddTrialIntegrator(mfem.TransposeIntegrator(mfem.MixedCurlIntegrator(one)),
                         None,
                         TrialSpace["E_space"],
                         TestSpace["F_space"])
    # -i ω ϵ (E , G) = i (- ω ϵ E, G)
    a.AddTrialIntegrator(None,
                         mfem.TransposeIntegrator(
                             mfem.VectorFEMassIntegrator(negepsomeg_cf)),
                         TrialSpace["E_space"],
                         TestSpace["G_space"])
    #  (H,∇ × G)
    a.AddTrialIntegrator(mfem.TransposeIntegrator(mfem.MixedCurlIntegrator(one)),
                         None,
                         TrialSpace["H_space"],
                         TestSpace["G_space"])
    # < n×Ĥ ,G>
    a.AddTrialIntegrator(mfem.TangentTraceIntegrator(), None,
                         TrialSpace["hatH_space"],
                         TestSpace["G_space"])
    # test integrators
    # (∇×G ,∇× δG)
    a.AddTestIntegrator(mfem.CurlCurlIntegrator(one), None,
                        TestSpace["G_space"],
                        TestSpace["G_space"])
    # (G,δG)
    a.AddTestIntegrator(mfem.VectorFEMassIntegrator(one), None,
                        TestSpace["G_space"],
                        TestSpace["G_space"])

    if dim == 3:
        # i ω μ (H, F)
        a.AddTrialIntegrator(None, mfem.TransposeIntegrator(
            mfem.VectorFEMassIntegrator(muomeg_cf)),
            TrialSpace["H_space"],
            TestSpace["F_space"])
        # < n×Ê,F>
        a.AddTrialIntegrator(mfem.TangentTraceIntegrator(), None,
                             TrialSpace["hatE_space"],
                             TestSpace["F_space"])

        # test integrators
        # (∇×F,∇×δF)
        a.AddTestIntegrator(mfem.CurlCurlIntegrator(one), None,
                            TestSpace["F_space"],
                            TestSpace["F_space"])
        # (F,δF)
        a.AddTestIntegrator(mfem.VectorFEMassIntegrator(one), None,
                            TestSpace["F_space"],
                            TestSpace["F_space"])
        # μ^2 ω^2 (F,δF)
        a.AddTestIntegrator(mfem.VectorFEMassIntegrator(mu2omeg2_cf), None,
                            TestSpace["F_space"],
                            TestSpace["F_space"])
        # -i ω μ (F,∇ × δG) = i (F, -ω μ ∇ × δ G)
        a.AddTestIntegrator(None, mfem.MixedVectorWeakCurlIntegrator(negmuomeg_cf),
                            TestSpace["F_space"],
                            TestSpace["G_space"])
        # -i ω ϵ (∇ × F, δG)
        a.AddTestIntegrator(None, mfem.MixedVectorCurlIntegrator(negepsomeg_cf),
                            TestSpace["F_space"],
                            TestSpace["G_space"])
        # i ω μ (∇ × G,δF)
        a.AddTestIntegrator(None, mfem.MixedVectorCurlIntegrator(muomeg_cf),
                            TestSpace["G_space"],
                            TestSpace["F_space"])
        # i ω ϵ (G, ∇ × δF )
        a.AddTestIntegrator(None, mfem.MixedVectorWeakCurlIntegrator(epsomeg_cf),
                            TestSpace["G_space"],
                            TestSpace["F_space"])
        # ϵ^2 ω^2 (G,δG)
        a.AddTestIntegrator(mfem.VectorFEMassIntegrator(eps2omeg2_cf), None,
                            TestSpace["G_space"],
                            TestSpace["G_space"])
    else:
        # i ω μ (H, F)
        a.AddTrialIntegrator(None, mfem.MixedScalarMassIntegrator(muomeg_cf),
                             TrialSpace["H_space"],
                             TestSpace["F_space"])
        # < n×Ê,F>
        a.AddTrialIntegrator(mfem.TraceIntegrator(), None,
                             TrialSpace["hatE_space"],
                             TestSpace["F_space"])
        # test integrators
        # (∇F,∇δF)
        a.AddTestIntegrator(mfem.DiffusionIntegrator(one), None,
                            TestSpace["F_space"],
                            TestSpace["F_space"])
        # (F,δF)
        a.AddTestIntegrator(mfem.MassIntegrator(one), None,
                            TestSpace["F_space"],
                            TestSpace["F_space"])
        # μ^2 ω^2 (F,δF)
        a.AddTestIntegrator(mfem.MassIntegrator(mu2omeg2_cf), None,
                            TestSpace["F_space"],
                            TestSpace["F_space"])
        # -i ω μ (F,∇ × δG) = i (F, -ω μ ∇ × δ G)
        a.AddTestIntegrator(None,
                            mfem.TransposeIntegrator(
                                mfem.MixedCurlIntegrator(negmuomeg_cf)),
                            TestSpace["F_space"],
                            TestSpace["G_space"])
        # -i ω ϵ (∇ × F, δG) = i (- ω ϵ A ∇ F,δG), A = [0 1; -1; 0]
        a.AddTestIntegrator(None, mfem.MixedVectorGradientIntegrator(negepsrot_cf),
                            TestSpace["F_space"],
                            TestSpace["G_space"])
        # i ω μ (∇ × G,δF) = i (ω μ ∇ × G, δF )
        a.AddTestIntegrator(None, mfem.MixedCurlIntegrator(muomeg_cf),
                            TestSpace["G_space"],
                            TestSpace["F_space"])
        # i ω ϵ (G, ∇ × δF ) =  i (ω ϵ G, A ∇ δF) = i ( G , ω ϵ A ∇ δF)
        a.AddTestIntegrator(None,
                            mfem.TransposeIntegrator(
                                mfem.MixedVectorGradientIntegrator(epsrot_cf)),
                            TestSpace["G_space"],
                            TestSpace["F_space"])
        # ϵ^2 ω^2 (G, δG)
        a.AddTestIntegrator(mfem.VectorFEMassIntegrator(eps2omeg2_cf), None,
                            TestSpace["G_space"],
                            TestSpace["G_space"])
    if with_pml:
        # trial integrators
        # -i ω ϵ (β E , G) = -i ω ϵ ((β_re + i β_im) E, G)
        #                  = (ω ϵ β_im E, G) + i (- ω ϵ β_re E, G)
        a.AddTrialIntegrator(mfem.TransposeIntegrator(mfem.VectorFEMassIntegrator(
            epsomeg_detJ_Jt_J_inv_i_restr)),
            mfem.TransposeIntegrator(mfem.VectorFEMassIntegrator(
                negepsomeg_detJ_Jt_J_inv_r_restr)),
            TrialSpace["E_space"], TestSpace["G_space"])
        if dim == 3:
            # trial integrators
            # i ω μ (α^-1 H, F) = i ω μ ( (α^-1_re + i α^-1_im) H, F)
            #                   = (- ω μ α^-1_im, H,F) + i *(ω μ α^-1_re H, F)
            a.AddTrialIntegrator(
                mfem.TransposeIntegrator(mfem.VectorFEMassIntegrator(
                    negmuomeg_detJ_Jt_J_inv_i_restr)),
                mfem.TransposeIntegrator(mfem.VectorFEMassIntegrator(
                    muomeg_detJ_Jt_J_inv_r_restr)),
                TrialSpace["H_space"], TestSpace["F_space"])
            # test integrators
            # μ^2 ω^2 (|α|^-2 F,δF)
            a.AddTestIntegrator(
                mfem.VectorFEMassIntegrator(
                    mu2omeg2_detJ_Jt_J_inv_2_restr), None,
                TestSpace["F_space"], TestSpace["F_space"])
            # -i ω μ (α^-* F,∇ × δG) = i (F, - ω μ α^-1 ∇ × δ G)
            #                        = i (F, - ω μ (α^-1_re + i α^-1_im) ∇ × δ G)
            #                        = (F, - ω μ α^-1_im ∇ × δ G) + i (F, - ω μ α^-1_re ∇×δG)
            a.AddTestIntegrator(mfem.MixedVectorWeakCurlIntegrator(
                negmuomeg_detJ_Jt_J_inv_i_restr),
                mfem.MixedVectorWeakCurlIntegrator(
                negmuomeg_detJ_Jt_J_inv_r_restr),
                TestSpace["F_space"], TestSpace["G_space"])
            # -i ω ϵ (β ∇ × F, δG) = -i ω ϵ ((β_re + i β_im) ∇ × F, δG)
            #                      = (ω ϵ β_im  ∇ × F, δG) + i (- ω ϵ β_re ∇ × F, δG)
            a.AddTestIntegrator(mfem.MixedVectorCurlIntegrator(
                epsomeg_detJ_Jt_J_inv_i_restr),
                mfem.MixedVectorCurlIntegrator(
                negepsomeg_detJ_Jt_J_inv_r_restr),
                TestSpace["F_space"], TestSpace["G_space"])
            # i ω μ (α^-1 ∇ × G,δF) = i ω μ ((α^-1_re + i α^-1_im) ∇ × G,δF)
            #                       = (- ω μ α^-1_im ∇ × G,δF) + i (ω μ α^-1_re ∇ × G,δF)
            a.AddTestIntegrator(mfem.MixedVectorCurlIntegrator(
                negmuomeg_detJ_Jt_J_inv_i_restr),
                mfem.MixedVectorCurlIntegrator(
                muomeg_detJ_Jt_J_inv_r_restr),
                TestSpace["G_space"], TestSpace["F_space"])
            # i ω ϵ (β^* G, ∇×δF) = i ω ϵ ( (β_re - i β_im) G, ∇×δF)
            #                     = (ω ϵ β_im G, ∇×δF) + i ( ω ϵ β_re G, ∇×δF)
            a.AddTestIntegrator(mfem.MixedVectorWeakCurlIntegrator(
                epsomeg_detJ_Jt_J_inv_i_restr),
                mfem.MixedVectorWeakCurlIntegrator(
                epsomeg_detJ_Jt_J_inv_r_restr),
                TestSpace["G_space"], TestSpace["F_space"])
            # ϵ^2 ω^2 (|β|^2 G,δG)
            a.AddTestIntegrator(mfem.VectorFEMassIntegrator(
                eps2omeg2_detJ_Jt_J_inv_2_restr), None,
                TestSpace["G_space"], TestSpace["G_space"])

        else:
            # trial integrators
            # i ω μ (α^-1 H, F) = i ω μ ( (α^-1_re + i α^-1_im) H, F)
            #                   = (- ω μ α^-1_im, H,F) + i *(ω μ α^-1_re H, F)
            a.AddTrialIntegrator(
                mfem.MixedScalarMassIntegrator(negmuomeg_detJ_i_restr),
                mfem.MixedScalarMassIntegrator(muomeg_detJ_r_restr),
                TrialSpace["H_space"], TestSpace["F_space"])
            # test integrators
            # μ^2 ω^2 (|α|^-2 F,δF)
            a.AddTestIntegrator(mfem.MassIntegrator(mu2omeg2_detJ_2_restr), None,
                                TestSpace["F_space"], TestSpace["F_space"])
            # -i ω μ (α^-* F,∇ × δG) = (F, ω μ α^-1 ∇ × δ G)
            #                        =(F, - ω μ α^-1_im ∇ × δ G) + i (F, - ω μ α^-1_re ∇×δG)
            a.AddTestIntegrator(
                mfem.TransposeIntegrator(
                    mfem.MixedCurlIntegrator(negmuomeg_detJ_i_restr)),
                mfem.TransposeIntegrator(
                    mfem.MixedCurlIntegrator(negmuomeg_detJ_r_restr)),
                TestSpace["F_space"], TestSpace["G_space"])
            # -i ω ϵ (β ∇ × F, δG) = i (- ω ϵ β A ∇ F,δG), A = [0 1; -1; 0]
            #                      = (ω ϵ β_im A ∇ F, δG) + i (- ω ϵ β_re A ∇ F, δG)
            a.AddTestIntegrator(mfem.MixedVectorGradientIntegrator(
                epsomeg_detJ_Jt_J_inv_i_rot_restr),
                mfem.MixedVectorGradientIntegrator(
                negepsomeg_detJ_Jt_J_inv_r_rot_restr),
                TestSpace["F_space"], TestSpace["G_space"])
            # i ω μ (α^-1 ∇ × G,δF) = i (ω μ α^-1 ∇ × G, δF )
            #                       = (- ω μ α^-1_im ∇ × G,δF) + i (ω μ α^-1_re ∇ × G,δF)
            a.AddTestIntegrator(mfem.MixedCurlIntegrator(negmuomeg_detJ_i_restr),
                                mfem.MixedCurlIntegrator(muomeg_detJ_r_restr),
                                TestSpace["G_space"], TestSpace["F_space"])
            # i ω ϵ (β^* G, ∇ × δF ) = i ( G , ω ϵ β A ∇ δF)
            #                        =  ( G , ω ϵ β_im A ∇ δF) + i ( G , ω ϵ β_re A ∇ δF)
            a.AddTestIntegrator(
                mfem.TransposeIntegrator(mfem.MixedVectorGradientIntegrator(
                    epsomeg_detJ_Jt_J_inv_i_rot_restr)),
                mfem.TransposeIntegrator(mfem.MixedVectorGradientIntegrator(
                    epsomeg_detJ_Jt_J_inv_r_rot_restr)),
                TestSpace["G_space"], TestSpace["F_space"])
            # ϵ^2 ω^2 (|β|^2 G,δG)
            a.AddTestIntegrator(mfem.VectorFEMassIntegrator(
                eps2omeg2_detJ_Jt_J_inv_2_restr), None,
                TestSpace["G_space"], TestSpace["G_space"])

    # RHS
    fncs = Functions(dim, dimc, omega, epsilon, mu, prob)
    if prob == 0:
        f_rhs_r = fncs.rhs_func_r
        f_rhs_i = fncs.rhs_func_i
        a.AddDomainLFIntegrator(mfem.VectorFEDomainLFIntegrator(f_rhs_r),
                                mfem.VectorFEDomainLFIntegrator(f_rhs_i),
                                TestSpace["G_space"])
    elif prob == 2:
        f_source = fncs.source_function
        a.AddDomainLFIntegrator(mfem.VectorFEDomainLFIntegrator(f_source),
                                None,
                                TestSpace["G_space"])

    hatEex_r = fncs.hatE_exact_r
    hatEex_i = fncs.hatE_exact_i

    if myid == 0:
        txt = "\n  Ref |" + "    Dofs    |" + "    ω    |"
        if exact_known:
            txt = txt + "  L2 Error  |" + "  Rate  |"

        txt = txt + "  Residual  |" + "  Rate  |" + " PCG it |"
        print(txt)

        if exact_known:
            print("-"*82)
        else:
            print("-"*60)

    res0 = 0.
    err0 = 0.
    dof0 = 0

    elements_to_refine = mfem.intArray()

    # visualizaiton socket
    E_out_r = None
    H_out_r = None

    E_r = mfem.ParGridFunction()
    E_i = mfem.ParGridFunction()
    H_r = mfem.ParGridFunction()
    H_i = mfem.ParGridFunction()

    if static_cond:
        a.EnableStaticCondensation()

    for it in range(pr+1):
        a.Assemble()

        ess_tdof_list = mfem.intArray()
        ess_bdr = mfem.intArray()

        if pmesh.bdr_attributes.Size() != 0:
            ess_bdr.SetSize(pmesh.bdr_attributes.Max())
            ess_bdr.Assign(1)
            hatE_fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)
            if with_pml:
                ess_bdr.Assign(0)
                ess_bdr[1] = 1

        # shift the ess_tdofs
        for j in range(ess_tdof_list.Size()):
            ess_tdof_list[j] += E_fes.GetTrueVSize() + H_fes.GetTrueVSize()

        offsets = mfem.intArray(5)
        offsets[0] = 0
        offsets[1] = E_fes.GetVSize()
        offsets[2] = H_fes.GetVSize()
        offsets[3] = hatE_fes.GetVSize()
        offsets[4] = hatH_fes.GetVSize()

        offsets.PartialSum()
        offsetsl = offsets.ToList()

        x = mfem.Vector([0.]*(2*offsetsl[-1]))
        if prob != 2:
            hatE_gf_r = mfem.ParGridFunction(hatE_fes, x, offsetsl[2])
            hatE_gf_i = mfem.ParGridFunction(
                hatE_fes, x, offsetsl[-1] + offsetsl[2])
            if dim == 3:
                hatE_gf_r.ProjectBdrCoefficientTangent(hatEex_r, ess_bdr)
                hatE_gf_i.ProjectBdrCoefficientTangent(hatEex_i, ess_bdr)
            else:
                hatE_gf_r.ProjectBdrCoefficientNormal(hatEex_r, ess_bdr)
                hatE_gf_i.ProjectBdrCoefficientNormal(hatEex_i, ess_bdr)
        Ah = mfem.OperatorPtr()
        X = mfem.Vector()
        B = mfem.Vector()

        a.FormLinearSystem(ess_tdof_list, x, Ah, X, B)

        Ahc = Ah.AsComplexOperator()

        BlockA_r = mfem.Opr2BlockOpr(Ahc.real())
        BlockA_i = mfem.Opr2BlockOpr(Ahc.imag())
        num_blocks = BlockA_r.NumRowBlocks()

        # this is to debug matrix
        # for i in range(num_blocks):
        #    for j in range(num_blocks):
        #       (mfem.Opr2HypreParMatrix(BlockA_r.GetBlock(i, j))).Print(str(i)+"_"+str(j)+"r")
        #       (mfem.Opr2HypreParMatrix(BlockA_i.GetBlock(i, j))).Print(str(i)+"_"+str(j)+"i")

        tdof_offsets = mfem.intArray(2*num_blocks+1)

        tdof_offsets[0] = 0
        skip = 0 if static_cond else 2
        k = 2 if static_cond else 0
        for i in range(num_blocks):
            tdof_offsets[i+1] = trial_fes[i+k].GetTrueVSize()
            tdof_offsets[num_blocks+i+1] = trial_fes[i+k].GetTrueVSize()
        tdof_offsets.PartialSum()

        blockA = mfem.BlockOperator(tdof_offsets)
        for i in range(num_blocks):
            for j in range(num_blocks):
                blockA.SetBlock(i, j, BlockA_r.GetBlock(i, j))
                blockA.SetBlock(i, j+num_blocks, BlockA_i.GetBlock(i, j), -1.0)
                blockA.SetBlock(i+num_blocks, j+num_blocks,
                                BlockA_r.GetBlock(i, j))
                blockA.SetBlock(i+num_blocks, j, BlockA_i.GetBlock(i, j))
        X.Assign(0.0)
        M = mfem.BlockDiagonalPreconditioner(tdof_offsets)

        if not static_cond:
            E_mat = mfem.Opr2HypreParMatrix(BlockA_r.GetBlock(0, 0))
            H_mat = mfem.Opr2HypreParMatrix(BlockA_r.GetBlock(1, 1))

            solver_E = mfem.HypreBoomerAMG(E_mat)
            solver_E.SetPrintLevel(0)
            solver_E.SetSystemsOptions(dim)
            solver_H = mfem.HypreBoomerAMG(H_mat)
            solver_H.SetPrintLevel(0)
            solver_H.SetSystemsOptions(dim)
            M.SetDiagonalBlock(0, solver_E)
            M.SetDiagonalBlock(1, solver_H)
            M.SetDiagonalBlock(num_blocks, solver_E)
            M.SetDiagonalBlock(num_blocks+1, solver_H)

        hatE_mat = mfem.Opr2HypreParMatrix(BlockA_r.GetBlock(skip, skip))
        solver_hatE = mfem.HypreAMS(hatE_mat, hatE_fes)
        solver_hatE.SetPrintLevel(0)

        hatH_mat = mfem.Opr2HypreParMatrix(BlockA_r.GetBlock(skip+1, skip+1))
        if dim == 2:
            solver_hatH = mfem.HypreBoomerAMG(hatH_mat)
            solver_hatH.SetPrintLevel(0)
        else:
            solver_hatH = mfem.HypreAMS(hatH_mat, hatH_fes)
            solver_hatH.SetPrintLevel(0)

        M.SetDiagonalBlock(skip, solver_hatE)
        M.SetDiagonalBlock(skip+1, solver_hatH)
        M.SetDiagonalBlock(skip+num_blocks, solver_hatE)
        M.SetDiagonalBlock(skip+num_blocks+1, solver_hatH)

        cg = mfem.CGSolver(MPI.COMM_WORLD)
        cg.SetRelTol(1e-6)
        cg.SetMaxIter(10000)
        cg.SetPrintLevel(0)
        cg.SetPreconditioner(M)
        cg.SetOperator(blockA)
        cg.Mult(B, X)

        # for i in range(num_blocks):
        #  delete &M.GetDiagonalBlock(i);

        num_iter = cg.GetNumIterations()

        a.RecoverFEMSolution(X, x)

        residuals = a.ComputeResidual(x)

        residual = residuals.Norml2()
        maxresidual = residuals.Max()
        globalresidual = residual * residual

        maxresidual = MPI.COMM_WORLD.allreduce(maxresidual, op=MPI.MAX)
        globalresidual = MPI.COMM_WORLD.allreduce(globalresidual, op=MPI.SUM)
        globalresidual = np.sqrt(globalresidual)

        E_r.MakeRef(E_fes, x, 0)
        E_i.MakeRef(E_fes, x, offsetsl[-1])

        H_r.MakeRef(H_fes, x, offsetsl[1])
        H_i.MakeRef(H_fes, x, offsetsl[-1]+offsetsl[1])

        dofs = 0
        for i in range(trial_fes.Size()):
            dofs += trial_fes[i].GlobalTrueVSize()

        if exact_known:
            E_ex_r = fncs.E_exact_r
            E_ex_i = fncs.E_exact_i
            H_ex_r = fncs.H_exact_r
            H_ex_i = fncs.H_exact_i
            E_err_r = E_r.ComputeL2Error(E_ex_r)
            E_err_i = E_r.ComputeL2Error(E_ex_i)
            H_err_r = H_r.ComputeL2Error(E_ex_r)
            H_err_i = H_r.ComputeL2Error(E_ex_i)
            L2Error = np.sqrt(E_err_r*E_err_r + E_err_i*E_err_i
                              + H_err_r*H_err_r + H_err_i*H_err_i)
            rate_err = 0 if it == 0 else dim * \
                np.log(err0/L2Error)/np.log(dof0/dofs)
            err0 = L2Error

        rate_res = 0.0 if it == 0 else dim * \
            np.log(res0/globalresidual)/np.log(dof0/dofs)

        res0 = globalresidual
        dof0 = dofs

        if myid == 0:
            txt = ("{:5d}".format(it) + " | " +
                   "{:10d}".format(dof0) + " | " +
                   " {:.1f}".format(2.0*rnum) + " π  | ")
            if exact_known:
                txt = txt + (" {:.3e}".format(err0) + " | " +
                             " {:.2f} ".format(rate_err) + " | ")

            txt = txt + (" {:.3e}".format(res0) + " | " +
                         " {:.2f} ".format(rate_res) + " | " +
                         "{:6d}".format(num_iter) + " | ")
            print(txt)

        if visualization:
            keys = "jRcml\n" if it == 0 and dim == 2 else ""
            E_out_r = VisualizeField(E_out_r, "localhost", visport, E_r,
                                     "Numerical Electric field (real part)", 0, 0, 500, 500, keys)
            H_out_r = VisualizeField(H_out_r, "localhost", visport, H_r,
                                     "Numerical Magnetic field (real part)", 0, 0, 500, 500, keys)

        if it == pr:
            break

        if theta > 0.0:
            elements_to_refine.SetSize(0)
            for iel in range(pmesh.GetNE()):
                if residuals[iel] > theta * maxresidual:
                    elements_to_refine.Append(iel)
            pmesh.GeneralRefinement(elements_to_refine, 1, 1)
        else:
            pmesh.UniformRefinement()

        if with_pml:
            pml.SetAttributes(pmesh)
        for i in range(trial_fes.Size()):
            trial_fes[i].Update(False)
        a.Update()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex40 (Eikonal queation)')
    parser.add_argument('-m', '--mesh',
                        default="inline-quad.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).")

    parser.add_argument("-rnum", "--number-of-wavelengths",
                        action='store', default=1.0, type=float,
                        help="Number of wavelengths")
    parser.add_argument("-mu", "--permeability",
                        action='store', default=1.0, type=float,
                        help="Permeability of free space (or 1/(spring constant)).")
    parser.add_argument("-eps", "--permittivity",
                        action='store', default=1.0, type=float,
                        help="Permittivity of free space (or mass constant).")
    parser.add_argument("-prob", "--problem",
                        action='store', default=0, type=int,
                        help="\n".join(("Problem case"
                                        " 0: plane wave, 1: Fichera 'oven', "
                                       " 2: Generic PML problem with point source given as a load "
                                        " 3: Scattering of a plane wave, "
                                        " 4: Point source given on the boundary")))
    parser.add_argument("-do", "--delta-order",
                        action='store', default=1, type=int,
                        help="Order enrichment for DPG test space.")
    parser.add_argument("-theta", "--theta",
                        action='store', default=0.0, type=float,
                        help="Theta parameter for AMR")
    parser.add_argument("-sref", "--serial-ref",
                        action='store', default=0, type=int,
                        help="Number of parallel refinements.")
    parser.add_argument("-pref", "--parallel-ref",
                        action='store', default=1, type=int,
                        help="Number of parallel refinements.")
    parser.add_argument("-sc", "--static-condensation",
                        action='store_true', default=False,
                        help="Enable static condensation.")
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        default=False,
                        help='Disable or disable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    visualization = not args.no_visualization

    run(meshfile=args.mesh,
        order=args.order,
        prob=args.problem,
        sr=args.serial_ref,
        pr=args.parallel_ref,
        epsilon=args.permittivity,
        mu=args.permeability,
        delta_order=args.delta_order,
        rnum=args.number_of_wavelengths,
        theta=args.theta,
        static_cond=args.static_condensation,
        visualization=visualization)
