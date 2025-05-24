##
##                   MFEM Ultraweak DPG Maxwell parallel example 
##
## Compile with: make pmaxwell
##
## sample run
##  mpirun -np 4 python pmaxwell.py -m ../../data/star.mesh -o 2 -sref 0 -pref 3 -rnum 0.5 -prob 0
##  mpirun -np 4 python pmaxwell.py -m ../../data/inline-quad.mesh -o 3 -sref 0 -pref 3 -rnum 4.8 -sc -prob 0
##  mpirun -np 4 python pmaxwell.py -m ../../data/inline-hex.mesh -o 2 -sref 0 -pref 1 -rnum 0.8 -sc -prob 0
##  mpirun -np 4 python pmaxwell.py -m ../../data/inline-quad.mesh -o 3 -sref 1 -pref 3 -rnum 4.8 -sc -prob 2
##  mpirun -np 4 python pmaxwell.py -o 3 -sref 1 -pref 2 -rnum 11.8 -sc -prob 3
##  mpirun -np 4 python pmaxwell.py -o 3 -sref 1 -pref 2 -rnum 9.8 -sc -prob 4

## AMR run. Note that this is a computationally intensive sample run.
## We recommend trying it on a large machine with more mpi ranks
##  mpirun -np 4 pmaxwell -o 3 -sref 0 -pref 15 -prob 1 -theta 0.7 -sc

## Description:
## This example code demonstrates the use of MFEM to define and solve
## the "ultraweak" (UW) DPG formulation for the Maxwell problem

##      ∇×(1/μ ∇×E) - ω² ϵ E = Ĵ ,   in Ω
##                       E×n = E₀ , on ∂Ω

## It solves the following kinds of problems
## 1) Known exact solutions with error convergence rates
##    a) A manufactured solution problem where E is a plane beam
## 2) Fichera "microwave" problem
## 3) PML problems
##    a) Generic PML problem with point source given by the load
##    b) Plane wave scattering from a square
##    c) PML problem with a point source prescribed on the boundary

## The DPG UW deals with the First Order System
##  i ω μ H + ∇ × E = 0,   in Ω
## -i ω ϵ E + ∇ × H = J,   in Ω
##            E × n = E_0, on ∂Ω
## Note: Ĵ = -iωJ

## The ultraweak-DPG formulation is obtained by integration by parts of both
## equations and the introduction of trace unknowns on the mesh skeleton

## in 2D
## E is vector valued and H is scalar.
##    (∇ × E, F) = (E, ∇ × F) + < n × E , F>
## or (∇ ⋅ AE , F) = (AE, ∇ F) + < AE ⋅ n, F>
## where A = [0 1; -1 0];

## E ∈ (L²(Ω))² , H ∈ L²(Ω)
## Ê ∈ H^-1/2(Γₕ), Ĥ ∈ H^1/2(Γₕ)
##  i ω μ (H,F) + (E, ∇ × F) + < AÊ, F > = 0,      ∀ F ∈ H¹
## -i ω ϵ (E,G) + (H,∇ × G)  + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
##                                        Ê = E₀      on ∂Ω
## -------------------------------------------------------------------------
## |   |       E      |      H      |      Ê       |       Ĥ      |  RHS    |
## -------------------------------------------------------------------------
## | F |  (E,∇ × F)   | i ω μ (H,F) |   < Ê, F >   |              |         |
## |   |              |             |              |              |         |
## | G | -i ω ϵ (E,G) |  (H,∇ × G)  |              | < Ĥ, G × n > |  (J,G)  |
## where (F,G) ∈  H¹ × H(curl,Ω)

## in 3D
## E,H ∈ (L^2(Ω))³
## Ê ∈ H_0^1/2(Ω)(curl, Γₕ), Ĥ ∈ H^-1/2(curl, Γₕ)
##  i ω μ (H,F) + (E,∇ × F) + < Ê, F × n > = 0,      ∀ F ∈ H(curl,Ω)
## -i ω ϵ (E,G) + (H,∇ × G) + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
##                                   Ê × n = E₀      on ∂Ω
## -------------------------------------------------------------------------
## |   |       E      |      H      |      Ê       |       Ĥ      |  RHS    |
## -------------------------------------------------------------------------
## | F |  (E,∇ × F)   | i ω μ (H,F) | < n × Ê, F > |              |         |
## |   |              |             |              |              |         |
## | G | -i ω ϵ (E,G) |  (H,∇ × G)  |              | < n × Ĥ, G > |  (J,G)  |
## where (F,G) ∈  H(curl,Ω) × H(curl,Ω)

## Here we use the "Adjoint Graph" norm on the test space i.e.,
## ||(F,G)||²ᵥ  = ||A^*(F,G)||² + ||(F,G)||² where A is the
## maxwell operator defined by (1)

## The PML formulation is

##      ∇×(1/μ α ∇×E) - ω² ϵ β E = Ĵ ,   in Ω
##                E×n = E₀ , on ∂Ω

## where α = |J|⁻¹ Jᵀ J (in 2D it's the scalar |J|⁻¹),
## β = |J| J⁻¹ J⁻ᵀ, J is the Jacobian of the stretching map
## and |J| its determinant.

## The first order system reads
##  i ω μ α⁻¹ H + ∇ × E = 0,   in Ω
##    -i ω ϵ β E + ∇ × H = J,   in Ω
##                 E × n = E₀,  on ∂Ω

## and the ultraweak formulation is

## in 2D
## E ∈ (L²(Ω))² , H ∈ L²(Ω)
## Ê ∈ H^-1/2(Ω)(Γₕ), Ĥ ∈ H^1/2(Γₕ)
##  i ω μ (α⁻¹ H,F) + (E, ∇ × F) + < AÊ, F > = 0,          ∀ F ∈ H¹
## -i ω ϵ (β E,G)   + (H,∇ × G)  + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
##                                            Ê = E₀     on ∂Ω
## ---------------------------------------------------------------------------------
## |   |       E        |        H         |      Ê       |       Ĥ      |  RHS    |
## ---------------------------------------------------------------------------------
## | F |  (E,∇ × F)     | i ω μ (α⁻¹ H,F)  |   < Ê, F >   |              |         |
## |   |                |                  |              |              |         |
## | G | -i ω ϵ (β E,G) |    (H,∇ × G)     |              | < Ĥ, G × n > |  (J,G)  |

## where (F,G) ∈  H¹ × H(curl,Ω)

##
## in 3D
## E,H ∈ (L^2(Ω))³
## Ê ∈ H_0^1/2(Ω)(curl, Γ_h), Ĥ ∈ H^-1/2(curl, Γₕ)
##  i ω μ (α⁻¹ H,F) + (E,∇ × F) + < Ê, F × n > = 0,      ∀ F ∈ H(curl,Ω)
## -i ω ϵ (β E,G)    + (H,∇ × G) + < Ĥ, G × n > = (J,G)   ∀ G ∈ H(curl,Ω)
##                                        Ê × n = E_0     on ∂Ω
## -------------------------------------------------------------------------------
## |   |       E      |      H           |      Ê       |       Ĥ      |  RHS    |
## -------------------------------------------------------------------------------
## | F |  ( E,∇ × F)  | i ω μ (α⁻¹ H,F)  | < n × Ê, F > |              |         |
## |   |              |                  |              |              |         |
## | G | -iωϵ (β E,G) |   (H,∇ × G)      |              | < n × Ĥ, G > |  (J,G)  |
## where (F,G) ∈  H(curl,Ω) × H(curl,Ω)

## For more information see https:##doi.org/10.1016/j.camwa.2021.01.017
import os
from os.path import expanduser, join

import numpy as np
from numpy import pi

import mfem.par as mfem

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)

prob_type = {"plane_wave",
             "fichera_oven",
             "pml_general",
             "pml_plane_wave_scatter",
             "pml_pointsource"}

def run(meshfile='',
        order=1,
        delta_order=1,        
        prob=0,
        sr=0, 
        pr=1, 
        eps=1.0,
        mu=1.0,
        rnum=1.0,
        theta=0.0,
        static_cond=False,
        visualization=False):

    omega = 2.*pi*rnum
    with_pml = False
    
    if prob == 0:
       exact_known = True;
    elif prob == 1:
       meshfile = "meshes/fichera-waveguide.mesh";
       omega = 5.0;
       rnum = omega/(2.*M_PI);
    elif prob == 2:
       with_pml = True
    else:
       with_pml = True
       meshfile = "meshes/scatter.mesh"
       
    mesh_file = expanduser(
        join(os.path.dirname(__file__), '..', 'data', meshfile))

    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()
    assert dim > 1, "Dimension = 1 is not supported in this example"

    dimc = 3 if dim == 3 else 1
    for i in range(sr):
        mesh.UniformRefinement();
    mesh.EnsureNCMesh(False)


    if with_pml:
        assert False, "PML is not supported"
        
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh

    attr = mfem.intArray()
    attrPML= mfem.intArray()    

    # Define spaces
    TrialSpace = {"E_space": 0,
                  "H_space": 1,
                  "hatE_space": 2,
                  "hatH_space": 3,}
    TestSpace = {"F_space": 0,
                 "G_space": 1,}

    # Vector L2 L2 space for E
    E_fec = mfem.FiniteElementCollection(order-1, dim)
    E_fes = mfem.ParFiniteElementSpace(pmesh, E_fec, dim)
    
    # Vector L2 L2 space for H
    H_fec = mfem.FiniteElementCollection(order-1, dim)
    H_fes = mfem.ParFiniteElementSpace(pmesh, H_fec, dim)

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
    


    trial_fes = mfem.ParFiniteElementSpacePtrArray()
    test_fec = mfem.FiniteElementCollectionPtrArray()
    trial_fes.Append(E_fes)
    trial_fes.Append(H_fes)
    trial_fes.Append(hatE_fes)
    trial_fes.Append(hatH_fes)
    test_fec.Append(F_fec)
    test_fec.Append(G_fec)
     
    # Bilinear form coefficients
    one = mfem.ConstantCoefficient(1.0);
    eps2omeg2 = mfem.ConstantCoefficient(epsilon*epsilon*omega*omega)
    mu2omeg2 = mfem.ConstantCoefficient(mu*mu*omega*omega)
    muomeg = mfem.ConstantCoefficient(mu*omega)
    negepsomeg = mfem.ConstantCoefficient(-epsilon*omega)
    epsomeg = mfem.ConstantCoefficient(epsilon*omega)
    negmuomeg = mfem.ConstantCoefficient(-mu*omega)

    # for the 2D case
    rot_mat = mfem.DenseMatrix(np.array([[0, 1.],[-1, 0]])
    rot =mfem.MatrixConstantCoefficient(rot_mat)
    epsrot = mfem.ScalarMatrixProductCoefficient(epsomeg,rot)
    negepsrot = ScalarMatrixProductCoefficient(negepsomeg,rot)

    if pml:
        epsomeg_cf = mfem.RestrictedCoefficient(epsomeg,attr)
        negepsomeg_cf = mfem.RestrictedCoefficient(negepsomeg,attr)
        eps2omeg2_cf = mfem.RestrictedCoefficient(eps2omeg2,attr)
        muomeg_cf = mfem.RestrictedCoefficient(muomeg,attr)
        negmuomeg_cf = mfem.RestrictedCoefficient(negmuomeg,attr)
        mu2omeg2_cf = mfem.RestrictedCoefficient(mu2omeg2,attr)
        epsrot_cf = mfem.MatrixRestrictedCoefficient(epsrot,attr)
        negepsrot_cf = mfem.MatrixRestrictedCoefficient(negepsrot,attr)
    else:
        epsomeg_cf = epsomeg
        negepsomeg_cf = negepsomeg
        eps2omeg2_cf = eps2omeg2
        muomeg_cf = muomeg
        negmuomeg_cf = negmuomeg
        mu2omeg2_cf = mu2omeg2
        epsrot_cf = epsrot
        negepsrot_cf = negepsrot
                               
    detJ_r = mfem.PmlCoefficient(detJ_r_function,pml)
    detJ_i = mfem.PmlCoefficient(detJ_i_function,pml)
    abs_detJ_2 = mfem.PmlCoefficient(abs_detJ_2_function,pml)
    detJ_Jt_J_inv_r = mfem.PmlMatrixCoefficient(dim,detJ_Jt_J_inv_r_function,pml)
    detJ_Jt_J_inv_i = mfem.PmlMatrixCoefficient(dim,detJ_Jt_J_inv_i_function,pml)
    abs_detJ_Jt_J_inv_2 = mfem.PmlMatrixCoefficient(dim,abs_detJ_Jt_J_inv_2_function,pml)
    negmuomeg_detJ_r = mfem.ProductCoefficient(negmuomeg,detJ_r)
    negmuomeg_detJ_i = mfem.ProductCoefficient(negmuomeg,detJ_i)
    muomeg_detJ_r = mfem.ProductCoefficient(muomeg,detJ_r)
    mu2omeg2_detJ_2 = mfem.ProductCoefficient(mu2omeg2,abs_detJ_2)
    epsomeg_detJ_Jt_J_inv_i = mfem.ScalarMatrixProductCoefficient(epsomeg, detJ_Jt_J_inv_i)
    epsomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(epsomeg, detJ_Jt_J_inv_r)
    negepsomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(negepsomeg, detJ_Jt_J_inv_r)
    muomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(muomeg,detJ_Jt_J_inv_r)
    negmuomeg_detJ_Jt_J_inv_i = mfem.ScalarMatrixProductCoefficient(negmuomeg, detJ_Jt_J_inv_i)
    negmuomeg_detJ_Jt_J_inv_r = mfem.ScalarMatrixProductCoefficient(negmuomeg, detJ_Jt_J_inv_r)
    mu2omeg2_detJ_Jt_J_inv_2 = mfem.ScalarMatrixProductCoefficient(mu2omeg2, abs_detJ_Jt_J_inv_2)
    eps2omeg2_detJ_Jt_J_inv_2 = mfem.ScalarMatrixProductCoefficient(eps2omeg2, abs_detJ_Jt_J_inv_2)
    negmuomeg_detJ_r_restr = mfem.RestrictedCoefficient(negmuomeg_detJ_r,attrPML)
    negmuomeg_detJ_i_restr = mfem.RestrictedCoefficient(negmuomeg_detJ_i,attrPML)
    muomeg_detJ_r_restr = mfem.RestrictedCoefficient(muomeg_detJ_r,attrPML)
    mu2omeg2_detJ_2_restr = mfem.RestrictedCoefficient(mu2omeg2_detJ_2,attrPML)
    epsomeg_detJ_Jt_J_inv_i_restr = mfem.MatrixRestrictedCoefficient(epsomeg_detJ_Jt_J_inv_i,attrPML)
    epsomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(epsomeg_detJ_Jt_J_inv_r,attrPML)
    negepsomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(negepsomeg_detJ_Jt_J_inv_r,attrPML)
    muomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(muomeg_detJ_Jt_J_inv_r, attrPML)
    negmuomeg_detJ_Jt_J_inv_i_restr = mfem.MatrixRestrictedCoefficient(negmuomeg_detJ_Jt_J_inv_i,attrPML)
    negmuomeg_detJ_Jt_J_inv_r_restr = mfem.MatrixRestrictedCoefficient(negmuomeg_detJ_Jt_J_inv_r,attrPML)
    mu2omeg2_detJ_Jt_J_inv_2_restr = mfem.MatrixRestrictedCoefficient(mu2omeg2_detJ_Jt_J_inv_2,attrPML)
    eps2omeg2_detJ_Jt_J_inv_2_restr = mfem.MatrixRestrictedCoefficient(eps2omeg2_detJ_Jt_J_inv_2,attrPML)
                               

    if pml and dim == 2:
        epsomeg_detJ_Jt_J_inv_i_rot = mfem.MatrixProductCoefficient(epsomeg_detJ_Jt_J_inv_i, rot)
        epsomeg_detJ_Jt_J_inv_r_rot = mfem.MatrixProductCoefficient(epsomeg_detJ_Jt_J_inv_r, rot)
        negepsomeg_detJ_Jt_J_inv_r_rot = mfem.MatrixProductCoefficient(negepsomeg_detJ_Jt_J_inv_r, rot)
        epsomeg_detJ_Jt_J_inv_i_rot_restr = mfem.MatrixRestrictedCoefficient(epsomeg_detJ_Jt_J_inv_i_rot,
                                                                             attrPML)
        epsomeg_detJ_Jt_J_inv_r_rot_restr = mfem.MatrixRestrictedCoefficient(epsomeg_detJ_Jt_J_inv_r_rot,
                                                                             attrPML)
        negepsomeg_detJ_Jt_J_inv_r_rot_restr = mfem.MatrixRestrictedCoefficient(negepsomeg_detJ_Jt_J_inv_r_rot,
                                                                                attrPML)


    a = mfem.dpg.ParComplexDPGWeakForm(trial_fes,test_fec)
    a.StoreMatrices()  # needed for AMR

    # (E,∇ × F)
    a.AddTrialIntegrator(mfem.TransposeIntegrator(mfem.MixedCurlIntegrator(one)),
                         None, 
                         TrialSpace["E_space"],
                         TestSpace["F_space"])
    # -i ω ϵ (E , G) = i (- ω ϵ E, G)
    a.AddTrialIntegrator(None,
                         mfem.TransposeIntegrator(mfem.VectorFEMassIntegrator(negepsomeg_cf)),
                         TrialSpace["E_space"],
                         TestSpace["G_space"])
    #  (H,∇ × G)
    a.AddTrialIntegrator(mfem.TransposeIntegrator(mfem.MixedCurlIntegrator(one)),
                         None,
                         TrialSpace["H_space"],
                         TestSpace["G_space"])
    # < n×Ĥ ,G>
    a.AddTrialIntegrator(mfem.TangentTraceIntegrator, None,
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
        a.AddTrialIntegrator(mfem.TangentTraceIntegrator, None,
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
        a.AddTrialIntegrator(mfem.TraceIntegrator, None,
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
                            mfem.TransposeIntegrator(mfem.MixedCurlIntegrator(negmuomeg_cf)),
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
                            mfem.TransposeIntegrator(mfem.MixedVectorGradientIntegrator(epsrot_cf))
                            TestSpace["G_space"],
                            TestSpace["F_space"])
        # ϵ^2 ω^2 (G, δG)
        a.AddTestIntegrator(mfem.VectorFEMassIntegrator(eps2omeg2_cf), None,
                            TestSpace["G_space"],
                            TestSpace["G_space"])
    if pml:
        assert False, "PML not supported yet"


    # RHS
    f_rhs_r = mfem.VectorFunctionCoefficient(dim, rhs_func_r)
    f_rhs_i = mfem.VectorFunctionCoefficient(dim, rhs_func_i)
    f_source = mfem.VectorFunctionCoefficient(dim, source_function)
                               
    if prob == 0:
        a.AddDomainLFIntegrator(mfem.VectorFEDomainLFIntegrator(f_rhs_r),
                                mfem.VectorFEDomainLFIntegrator(f_rhs_i),
                                TestSpace["G_space"])
    elif prob == 2:
        a.AddDomainLFIntegrator(mfem.VectorFEDomainLFIntegrator(f_source),
                                None,
                                TestSpace["G_space"])
 
    hatEex_r = mfem.VectorFunctionCoefficient(dim, hatE_exact_r)
    hatEex_i = mfem.VectorFunctionCoefficient(dim, hatE_exact_i)

    if myid == 0:
       txt = "\n  Ref |" +  "    Dofs    |" +  "    ω    |"
       if exact_known:
          txt = txt + "  L2 Error  |" + "  Rate  |"
                               
       txt = txt + "  Residual  |" +  "  Rate  |" +  " PCG it |"
       print(txt)
                               
       if exact_known:
           print("-"*82)
       else:
           print("-"*60)
                               
    res0 = 0.
    err0 = 0.
    dof0 = 0
                               
    elements_to_refine = mfem.intArray()
                               
   socketstream E_out_r;
   socketstream H_out_r;


   Array<int> 

   ParGridFunction E_r, E_i, H_r, H_i;

   ParaViewDataCollection * paraview_dc = nullptr;

    if static_cond:
        a.EnableStaticCondensation()
   for (int it = 0; it<=pr; it++)
   {
      a->Assemble();

      Array<int> ess_tdof_list;
      Array<int> ess_bdr;
      if (pmesh.bdr_attributes.Size())
      {
         ess_bdr.SetSize(pmesh.bdr_attributes.Max());
         ess_bdr = 1;
         hatE_fes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
         if (pml)
         {
            ess_bdr = 0;
            ess_bdr[1] = 1;
         }
      }

      // shift the ess_tdofs
      for (int j = 0; j < ess_tdof_list.Size(); j++)
      {
         ess_tdof_list[j] += E_fes->GetTrueVSize() + H_fes->GetTrueVSize();
      }

      Array<int> offsets(5);
      offsets[0] = 0;
      offsets[1] = E_fes->GetVSize();
      offsets[2] = H_fes->GetVSize();
      offsets[3] = hatE_fes->GetVSize();
      offsets[4] = hatH_fes->GetVSize();
      offsets.PartialSum();

      Vector x(2*offsets.Last());
      x = 0.;

      if (prob != 2)
      {
         ParGridFunction hatE_gf_r(hatE_fes, x, offsets[2]);
         ParGridFunction hatE_gf_i(hatE_fes, x, offsets.Last() + offsets[2]);
         if (dim == 3)
         {
            hatE_gf_r.ProjectBdrCoefficientTangent(hatEex_r, ess_bdr);
            hatE_gf_i.ProjectBdrCoefficientTangent(hatEex_i, ess_bdr);
         }
         else
         {
            hatE_gf_r.ProjectBdrCoefficientNormal(hatEex_r, ess_bdr);
            hatE_gf_i.ProjectBdrCoefficientNormal(hatEex_i, ess_bdr);
         }
      }

      OperatorPtr Ah;
      Vector X,B;
      a->FormLinearSystem(ess_tdof_list,x,Ah, X,B);

      ComplexOperator * Ahc = Ah.As<ComplexOperator>();

      BlockOperator * BlockA_r = dynamic_cast<BlockOperator *>(&Ahc->real());
      BlockOperator * BlockA_i = dynamic_cast<BlockOperator *>(&Ahc->imag());

      int num_blocks = BlockA_r->NumRowBlocks();
      Array<int> tdof_offsets(2*num_blocks+1);

      tdof_offsets[0] = 0;
      int skip = (static_cond) ? 0 : 2;
      int k = (static_cond) ? 2 : 0;
      for (int i=0; i<num_blocks; i++)
      {
         tdof_offsets[i+1] = trial_fes[i+k]->GetTrueVSize();
         tdof_offsets[num_blocks+i+1] = trial_fes[i+k]->GetTrueVSize();
      }
      tdof_offsets.PartialSum();

      BlockOperator blockA(tdof_offsets);
      for (int i = 0; i<num_blocks; i++)
      {
         for (int j = 0; j<num_blocks; j++)
         {
            blockA.SetBlock(i,j,&BlockA_r->GetBlock(i,j));
            blockA.SetBlock(i,j+num_blocks,&BlockA_i->GetBlock(i,j), -1.0);
            blockA.SetBlock(i+num_blocks,j+num_blocks,&BlockA_r->GetBlock(i,j));
            blockA.SetBlock(i+num_blocks,j,&BlockA_i->GetBlock(i,j));
         }
      }

      X = 0.;
      BlockDiagonalPreconditioner M(tdof_offsets);

      if (!static_cond)
      {
         HypreBoomerAMG * solver_E = new HypreBoomerAMG((HypreParMatrix &)
                                                        BlockA_r->GetBlock(0,0));
         solver_E->SetPrintLevel(0);
         solver_E->SetSystemsOptions(dim);
         HypreBoomerAMG * solver_H = new HypreBoomerAMG((HypreParMatrix &)
                                                        BlockA_r->GetBlock(1,1));
         solver_H->SetPrintLevel(0);
         solver_H->SetSystemsOptions(dim);
         M.SetDiagonalBlock(0,solver_E);
         M.SetDiagonalBlock(1,solver_H);
         M.SetDiagonalBlock(num_blocks,solver_E);
         M.SetDiagonalBlock(num_blocks+1,solver_H);
      }

      HypreSolver * solver_hatH = nullptr;
      HypreAMS * solver_hatE = new HypreAMS((HypreParMatrix &)BlockA_r->GetBlock(skip,
                                                                                 skip),
                                            hatE_fes);
      solver_hatE->SetPrintLevel(0);
      if (dim == 2)
      {
         solver_hatH = new HypreBoomerAMG((HypreParMatrix &)BlockA_r->GetBlock(skip+1,
                                                                               skip+1));
         dynamic_cast<HypreBoomerAMG*>(solver_hatH)->SetPrintLevel(0);
      }
      else
      {
         solver_hatH = new HypreAMS((HypreParMatrix &)BlockA_r->GetBlock(skip+1,skip+1),
                                    hatH_fes);
         dynamic_cast<HypreAMS*>(solver_hatH)->SetPrintLevel(0);
      }

      M.SetDiagonalBlock(skip,solver_hatE);
      M.SetDiagonalBlock(skip+1,solver_hatH);
      M.SetDiagonalBlock(skip+num_blocks,solver_hatE);
      M.SetDiagonalBlock(skip+num_blocks+1,solver_hatH);

      CGSolver cg(MPI_COMM_WORLD);
      cg.SetRelTol(1e-6);
      cg.SetMaxIter(10000);
      cg.SetPrintLevel(0);
      cg.SetPreconditioner(M);
      cg.SetOperator(blockA);
      cg.Mult(B, X);

      for (int i = 0; i<num_blocks; i++)
      {
         delete &M.GetDiagonalBlock(i);
      }

      int num_iter = cg.GetNumIterations();

      a->RecoverFEMSolution(X,x);

      Vector & residuals = a->ComputeResidual(x);

      real_t residual = residuals.Norml2();
      real_t maxresidual = residuals.Max();
      real_t globalresidual = residual * residual;
      MPI_Allreduce(MPI_IN_PLACE, &maxresidual, 1, MPITypeMap<real_t>::mpi_type,
                    MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &globalresidual, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_SUM, MPI_COMM_WORLD);

      globalresidual = sqrt(globalresidual);

      E_r.MakeRef(E_fes,x, 0);
      E_i.MakeRef(E_fes,x, offsets.Last());

      H_r.MakeRef(H_fes,x, offsets[1]);
      H_i.MakeRef(H_fes,x, offsets.Last()+offsets[1]);

      int dofs = 0;
      for (int i = 0; i<trial_fes.Size(); i++)
      {
         dofs += trial_fes[i]->GlobalTrueVSize();
      }

      real_t L2Error = 0.0;
      real_t rate_err = 0.0;

      if (exact_known)
      {
         VectorFunctionCoefficient E_ex_r(dim,E_exact_r);
         VectorFunctionCoefficient E_ex_i(dim,E_exact_i);
         VectorFunctionCoefficient H_ex_r(dim,H_exact_r);
         VectorFunctionCoefficient H_ex_i(dim,H_exact_i);
         real_t E_err_r = E_r.ComputeL2Error(E_ex_r);
         real_t E_err_i = E_i.ComputeL2Error(E_ex_i);
         real_t H_err_r = H_r.ComputeL2Error(H_ex_r);
         real_t H_err_i = H_i.ComputeL2Error(H_ex_i);
         L2Error = sqrt(  E_err_r*E_err_r + E_err_i*E_err_i
                          + H_err_r*H_err_r + H_err_i*H_err_i );
         rate_err = (it) ? dim*log(err0/L2Error)/log((real_t)dof0/dofs) : 0.0;
         err0 = L2Error;
      }

      real_t rate_res = (it) ? dim*log(res0/globalresidual)/log((
                                                                   real_t)dof0/dofs) : 0.0;

      res0 = globalresidual;
      dof0 = dofs;

      if (myid == 0)
      {
         std::ios oldState(nullptr);
         oldState.copyfmt(std::cout);
         std::cout << std::right << std::setw(5) << it << " | "
                   << std::setw(10) <<  dof0 << " | "
                   << std::setprecision(1) << std::fixed
                   << std::setw(4) <<  2.0*rnum << " π  | "
                   << std::setprecision(3);
         if (exact_known)
         {
            std::cout << std::setw(10) << std::scientific <<  err0 << " | "
                      << std::setprecision(2)
                      << std::setw(6) << std::fixed << rate_err << " | " ;
         }
         std::cout << std::setprecision(3)
                   << std::setw(10) << std::scientific <<  res0 << " | "
                   << std::setprecision(2)
                   << std::setw(6) << std::fixed << rate_res << " | "
                   << std::setw(6) << std::fixed << num_iter << " | "
                   << std::endl;
         std::cout.copyfmt(oldState);
      }

      if (visualization)
      {
         const char * keys = (it == 0 && dim == 2) ? "jRcml\n" : nullptr;
         char vishost[] = "localhost";
         VisualizeField(E_out_r,vishost, visport, E_r,
                        "Numerical Electric field (real part)", 0, 0, 500, 500, keys);
         VisualizeField(H_out_r,vishost, visport, H_r,
                        "Numerical Magnetic field (real part)", 501, 0, 500, 500, keys);
      }

      if (paraview)
      {
         paraview_dc->SetCycle(it);
         paraview_dc->SetTime((real_t)it);
         paraview_dc->Save();
      }

      if (it == pr)
      {
         break;
      }

      if (theta > 0.0)
      {
         elements_to_refine.SetSize(0);
         for (int iel = 0; iel<pmesh.GetNE(); iel++)
         {
            if (residuals[iel] > theta * maxresidual)
            {
               elements_to_refine.Append(iel);
            }
         }
         pmesh.GeneralRefinement(elements_to_refine,1,1);
      }
      else
      {
         pmesh.UniformRefinement();
      }
      if (pml) { pml->SetAttributes(&pmesh); }
      for (int i =0; i<trial_fes.Size(); i++)
      {
         trial_fes[i]->Update(false);
      }
      a->Update();
   }

                               
if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex40 (Eikonal queation)')
    parser.add_argument('-m', '--mesh',
                        default='star.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree).")
    
    parser.add_argument("-rnum", "--number-of-wavelengths",
                        action='store', default=1.0, type=float,
                        help="Number of wavelengths");
    parser.add_argument("-mu", "--permeability",
                        action='store', default=1.0, type=float, 
                        help="Permeability of free space (or 1/(spring constant)).")
    parser.add_argument("-eps", "--permittivity",
                        action='store', default=1.0, type=float,                         
                        help="Permittivity of free space (or mass constant).");
    parser.add_argument("-prob", "--problem",
                       action='store', default=1, type=int,
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
        eps=args.permittivity,
        mu=args.permeability,
        delta_order=args.delta_order,
        rnum=args.number_of_wavelengths,
        theta=args.theta,
        static_cond=args.static_condensation,
        visualization=visualization)

