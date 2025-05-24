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
                               

    pass

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

