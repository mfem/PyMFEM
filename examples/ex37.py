'''
   PyMFEM example 37

   See c++ version in the MFEM library for more detail

   Sample runs:
       python ex37.py -alpha 10
       python ex37.py -alpha 10 -pv
       python ex37.py -lambda 0.1 -mu 0.1
       python ex37.py -o 2 -alpha 5.0 -mi 50 -vf 0.4 -ntol 1e-5
       python ex37.py -r 6 -o 1 -alpha 25.0 -epsilon 0.02 -mi 50 -ntol 1e-5


   Description: This example code demonstrates the use of MFEM to solve a
                density-filtered [3] topology optimization problem. The
                objective is to minimize the compliance

                    minimize ∫_Ω f⋅u dx over u ∈ [H¹(Ω)]² and ρ ∈ L¹(Ω)

                    subject to

                      -Div(r(ρ̃)Cε(u)) = f       in Ω + BCs
                      -ϵ²Δρ̃ + ρ̃ = ρ             in Ω + Neumann BCs
                      0 ≤ ρ ≤ 1                 in Ω
                      ∫_Ω ρ dx = θ vol(Ω)

                Here, r(ρ̃) = ρ₀ + ρ̃³ (1-ρ₀) is the solid isotropic material
                penalization (SIMP) law, C is the elasticity tensor for an
                isotropic linearly elastic material, ϵ > 0 is the design
                length scale, and 0 < θ < 1 is the volume fraction.

                The problem is discretized and gradients are computing using
                finite elements [1]. The design is optimized using an entropic
                mirror descent algorithm introduced by Keith and Surowiec [2]
                that is tailored to the bound constraint 0 ≤ ρ ≤ 1.

                This example highlights the ability of MFEM to deliver high-
                order solutions to inverse design problems and showcases how
                to set up and solve PDE-constrained optimization problems
                using the so-called reduced space approach.

   [1] Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O.
       (2011). Efficient topology optimization in MATLAB using 88 lines of
       code. Structural and Multidisciplinary Optimization, 43(1), 1-16.
   [2] Keith, B. and Surowiec, T. (2023) Proximal Galerkin: A structure-
       preserving finite element method for pointwise bound constraints.
       arXiv:2307.12444 [math.NA]
   [3] Lazarov, B. S., & Sigmund, O. (2011). Filters in topology optimization
       based on Helmholtz‐type differential equations. International Journal
       for Numerical Methods in Engineering, 86(6), 765-781.

'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray

from ex37_common import (sigmoid,
                         der_sigmoid,
                         inv_sigmoid,
                         MappedGridFunctionCoefficient,
                         DiffMappedGridFunctionCoefficient,
                         SIMPInterpolationCoefficient,
                         StrainEnergyDensityCoefficient,
                         VolumeForceCoefficient,
                         DiffusionSolver,
                         LinearElasticitySolver,
                         )
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor

visualization = True
paraview_output = False


def proj(psi, target_volume, tol=1e-12, max_its=10):
    '''

    @brief Bregman projection of ρ = sigmoid(ψ) onto the subspace
          ∫_Ω ρ dx = θ vol(Ω) as follows:

          1. Compute the root of the R → R function
              f(c) = ∫_Ω sigmoid(ψ + c) dx - θ vol(Ω)
          2. Set ψ ← ψ + c.

    @param psi a GridFunction to be updated
    @param target_volume θ vol(Ω)
    @param tol Newton iteration tolerance
    @param max_its Newton maximum iteration number
    @return double Final volume, ∫_Ω sigmoid(ψ)

    '''

    sigmoid_psi = MappedGridFunctionCoefficient(psi, sigmoid)
    der_sigmoid_psi = MappedGridFunctionCoefficient(psi, der_sigmoid)

    int_sigmoid_psi = mfem.LinearForm(psi.FESpace())
    int_sigmoid_psi.AddDomainIntegrator(mfem.DomainLFIntegrator(sigmoid_psi))
    int_der_sigmoid_psi = mfem.LinearForm(psi.FESpace())
    int_der_sigmoid_psi.AddDomainIntegrator(mfem.DomainLFIntegrator(
                                              der_sigmoid_psi))
    done = False
    for k in range(max_its):  # Newton iteration
        int_sigmoid_psi.Assemble()         # Recompute f(c) with updated ψ
        f=int_sigmoid_psi.Sum() - target_volume

        int_der_sigmoid_psi.Assemble()      # Recompute df(c) with updated ψ
        df=int_der_sigmoid_psi.Sum()

        dc=-f/df
        psi += dc
        if abs(dc) < tol:
            done=True
            break

    if not done:
        message=("Projection reached maximum iteration without converging. " +
                   "Result may not be accurate.")
        import warnigns
        warnings.warn(message, RuntimeWarning)

    int_sigmoid_psi.Assemble();
    return int_sigmoid_psi.Sum();


def run(ref_levels = 5,
        order=2,
        alpha=1.0,
        epsilon=0.01,
        vol_fraction=0.5,
        max_it=1e3,
        itol=1e-1,
        ntol=1e-4,
        rho_min=1e-6,
        llambda=1.0,
        mu=1.0,):
    
    mesh = mfem.Mesh.MakeCartesian2D(3, 1, mfem.Element.QUADRILATERAL,
                                    True, 3.0, 1.0)
    dim = mesh.Dimension()

    # 2. Set BCs.
    for i in range(mesh.GetNBE()):
        be = mesh.GetBdrElement(i)
        vertices = mesh.GetBdrElementVertices(i)   # this method returns list

        coords1 = mesh.GetVertexArray(vertices[0])   # this returns numpy array
        coords2 = mesh.GetVertexArray(vertices[1])

        center = (coords1 + coords2)/2

        if abs(center[0] - 0.0) < 1e-10:
            # the left edge
            be.SetAttribute(1)
        else:
            # all other boundaries
            be.SetAttribute(2)
    mesh.SetAttributes();

    # 3. Refine the mesh.
    for lev in range(ref_levels):
       mesh.UniformRefinement()


    # 4. Define the necessary finite element spaces on the mesh.
    state_fec = mfem.H1_FECollection(order, dim)    # space for u
    filter_fec = mfem.H1_FECollection(order, dim)   # space for ρ̃
    control_fec = mfem.L2_FECollection (order-1, dim,
                                        mfem.BasisType.GaussLobatto)  # space for ψ
    state_fes = mfem.FiniteElementSpace(mesh, state_fec,dim)
    filter_fes = mfem.FiniteElementSpace(mesh, filter_fec)
    control_fes = mfem.FiniteElementSpace(mesh, control_fec)

    state_size = state_fes.GetTrueVSize()
    control_size = control_fes.GetTrueVSize()
    filter_size = filter_fes.GetTrueVSize()

    print("Number of state unknowns: " + str(state_size))
    print("Number of filter unknowns: " + str(filter_size))
    print("Number of control unknowns: " + str(control_size))
       
    #  5. Set the initial guess for ρ.
    u = mfem.GridFunction(state_fes)
    psi = mfem.GridFunction(control_fes)
    psi_old = mfem.GridFunction(control_fes)
    rho_filter = mfem.GridFunction(filter_fes)
    
    u.Assign(0.0)
    rho_filter.Assign(vol_fraction)
    psi.Assign(inv_sigmoid(vol_fraction))
    psi_old.Assign(inv_sigmoid(vol_fraction))

    # ρ = sigmoid(ψ)
    rho = MappedGridFunctionCoefficient(psi, sigmoid)
    #  Interpolation of ρ = sigmoid(ψ) in control fes (for ParaView output)
    rho_gf = mfem.GridFunction(control_fes)
    # ρ - ρ_old = sigmoid(ψ) - sigmoid(ψ_old)
    succ_diff_rho = DiffMappedGridFunctionCoefficient(psi, psi_old, sigmoid)

    # 6. Set-up the physics solver.
    maxat = mesh.bdr_attributes.Max();
    ess_bdr = intArray([0]*maxat)
    ess_bdr[0] = 1

    one = mfem.ConstantCoefficient(1.0)
    lambda_cf = mfem.ConstantCoefficient(llambda)
    mu_cf = mfem.ConstantCoefficient(mu)
    
    ElasticitySolver = LinearElasticitySolver()
    ElasticitySolver.mesh = mesh
    ElasticitySolver.order = state_fec.GetOrder()
    ElasticitySolver.SetupFEM()

    center = np.array((2.9, 0.5))
    force = np.array((0.0, -1.0))
    r = 0.05
    vforce_cf = VolumeForceCoefficient(r, center, force)
    ElasticitySolver.rhs_cf = vforce_cf
    ElasticitySolver.ess_bdr = ess_bdr

    # 7. Set-up the filter solver.
    eps2_cf = mfem.ConstantCoefficient(epsilon*epsilon)
    FilterSolver = DiffusionSolver()
    FilterSolver.mesh = mesh
    FilterSolver.order = filter_fec.GetOrder()
    FilterSolver.diffcf = eps2_cf
    FilterSolver.masscf = one
    
    if mesh.bdr_attributes.Size() > 0:
        ess_bdr_filter = mfem.intArray([0]*mesh.bdr_attributes.Max())
    else:
        ess_bdr_filter = mfem.intArray()

    FilterSolver.ess_bdr = ess_bdr_filter
    FilterSolver.SetupFEM()

    mass = mfem.BilinearForm(control_fes)
    mass.AddDomainIntegrator(mfem.InverseIntegrator(mfem.MassIntegrator(one)))
    mass.Assemble()
    M = mfem.SparseMatrix()
    empty = mfem.intArray()
    mass.FormSystemMatrix(empty,M)

    # 8. Define the Lagrange multiplier and gradient functions.
    grad = mfem.GridFunction(control_fes)
    w_filter = mfem.GridFunction(filter_fes)

    # 9. Define some tools for later.
    zero = mfem.ConstantCoefficient(0.0)
    onegf = mfem.GridFunction(control_fes)
    onegf.Assign(1.0)
    zerogf = mfem.GridFunction(control_fes)
    zerogf.Assign(0.0)
    
    vol_form = mfem.LinearForm(control_fes)
    vol_form.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
    vol_form.Assemble()
    domain_volume = vol_form(onegf)
    target_volume = domain_volume * vol_fraction

    # 10. Connect to GLVis. Prepare for VisIt output.
    if visualization:
        sout_r = mfem.socketstream("localhost", 19916) 
        sout_r.precision(8)

    if paraview_output:
       paraview_dc = mfem.ParaViewDataCollection("ex37", mesh)                               

       rho_gf.ProjectCoefficient(rho)
       paraview_dc.SetPrefixPath("ParaView")
       paraview_dc.SetLevelsOfDetail(order);
       paraview_dc.SetDataFormat(mfem.VTKFormat.BINARY)
       paraview_dc.SetHighOrderOutput(True)
       paraview_dc.SetCycle(0)
       paraview_dc.SetTime(0.0)
       paraview_dc.RegisterField("displacement", u)
       paraview_dc.RegisterField("density", rho_gf)
       paraview_dc.RegisterField("filtered_density", rho_filter)
       paraview_dc.Save();
 
           
    # 11. Iterate:
    for k in range(max_it):
        if k > 1:
            alpha *= k/(k-1.)

        print("\nStep = " + str(k))

        # Step 1 - Filter solve
        # Solve (ϵ^2 ∇ ρ̃, ∇ v ) + (ρ̃,v) = (ρ,v)
        FilterSolver.rhs_cf = rho
        FilterSolver.Solve();
        rho_filter = FilterSolver.GetFEMSolution()

        # Step 2 - State solve
        # Solve (λ r(ρ̃) ∇⋅u, ∇⋅v) + (2 μ r(ρ̃) ε(u), ε(v)) = (f,v)
        SIMP_cf = SIMPInterpolationCoefficient(rho_filter, rho_min, 1.0)
        lambda_SIMP_cf = mfem.ProductCoefficient(lambda_cf,SIMP_cf)
        mu_SIMP_cf = mfem.ProductCoefficient(mu_cf,SIMP_cf)
        
        ElasticitySolver.lambda_cf = lambda_SIMP_cf
        ElasticitySolver.mu_cf = mu_SIMP_cf
        
        ElasticitySolver.Solve()
        u = ElasticitySolver.GetFEMSolution()

        # Step 3 - Adjoint filter solve
        # Solve (ϵ² ∇ w̃, ∇ v) + (w̃ ,v) = (-r'(ρ̃) ( λ |∇⋅u|² + 2 μ |ε(u)|²),v)
        rhs_cf = StrainEnergyDensityCoefficient(lambda_cf, mu_cf, u, rho_filter,
                                                rho_min)
        FilterSolver.rhs_cf
        FilterSolver.Solve()
        w_filter = FilterSolver.GetFEMSolution()

        # Step 4 - Compute gradient
        # Solve G = M⁻¹w̃
        w_cf = mfem.GridFunctionCoefficient(w_filter)
        w_rhs = mfem.LinearForm(control_fes)
        w_rhs.AddDomainIntegrator(mfem.DomainLFIntegrator(w_cf))
        w_rhs.Assemble()
        M.Mult(w_rhs, grad)

        # Step 5 - Update design variable ψ ← proj(ψ - αG)
        psi.Add(-alpha, grad)
        material_volume = proj(psi, target_volume)

        # Compute ||ρ - ρ_old|| in control fes.
        norm_increment = zerogf.ComputeL1Error(succ_diff_rho)
        norm_reduced_gradient = norm_increment/alpha;
        psi_old.Assign(psi)

        compliance = (ElasticitySolver.GetLinearForm())(u)
        print("norm of the reduced gradient = " + "{g:}".format(norm_reduced_gradient))
        print("norm of the increment = " + "{g:}".format(norm_increment))
        print("compliance = " + "{g:}".format(compliance))
        print("volume fraction = " + "{g:}".format(material_volume / domain_volume))

        if visualization:
            r_gf = mfem.GridFunction(filter_fes)
            r_gf.ProjectCoefficient(SIMP_cf);
            sout_r << "solution\n" << mesh << r_gf << "window_title 'Design density r(ρ̃)'"
            sout_r.flush()
                
        if paraview_output:
           rho_gf.ProjectCoefficient(rho)
           paraview_dc.SetCycle(k)
           paraview_dc.SetTime(float(k))
           paraview_dc.Save()

        if norm_reduced_gradient < ntol and norm_increment < itol:
           break;

                      
if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex37 (Topology Optimization)')

    parser.add_argument('-r', '--ref_levels',
                        action='store', default=1, type=int,
                        help="Number of times to refine the mesh uniformly.")
    parser.add_argument("-o", "--order",
                        action='store', default=1, type=int,
                        help="Order (degree) of the finite elements.")
    parser.add_argument("-alpha", "--alpha-step-length",
                        action='store', default=1.0, type=float,
                        help="Step length for gradient descent.")
    parser.add_argument("-epsilon", "--epsilon-thickness",
                        action='store', default=0.01, type=float,
                        help="Length scale for ρ.")
    parser.add_argument("-mi", "--max-it",
                        action='store', default=1000, type=int,
                        help="Maximum number of gradient descent iterations.");
    parser.add_argument("-ntol", "--rel-tol",
                        action='store', default=1e-4, type=float,
                        help="Normalized exit tolerance.")
    parser.add_argument("-itol", "--abs-tol",
                        action='store', default=1e-1, type=float,
                        help="Increment exit tolerance.")
    parser.add_argument("-vf", "--volume-fraction",
                        action='store', default=0.5, type=float,
                        help="Volume fraction for the material density.")
    parser.add_argument("-lambda", "--llambda",
                        action='store', default=1.0, type=float,
                        help="Lamé constant λ.")
    parser.add_argument("-mu", "--mu",
                        action='store', default=1.0, type=float,
                        help="Lamé constant μ.")
    parser.add_argument("-rmin", "--psi-min",
                        action='store', default=1e-6, type=float,
                        help="Minimum of density coefficient.")
    parser.add_argument("-pv", "--paraview",
                        action='store_true', default=False,
                        help="Enable or disable ParaView output.");
    parser.add_argument("-no-vis", "--no-visualization",
                        action='store_false', default=True,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    globals()["visualization"] = args.no_visualization
    globals()["paraview_output"] = args.paraview

    run(ref_levels=args.ref_levels,
        order=args.order,
        alpha=args.alpha_step_length,
        epsilon=args.epsilon_thickness,
        vol_fraction=args.volume_fraction,
        max_it=args.max_it,
        itol=args.abs_tol,
        ntol=args.rel_tol,
        rho_min=args.psi_min,
        llambda=args.llambda,
        mu=args.mu)

