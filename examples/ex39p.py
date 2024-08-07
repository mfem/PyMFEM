'''
   MFEM example 39 (converted from ex39.cpp)

   See c++ version in the MFEM library for more detail

   Sample runs: mpirun -np 4 python ex39p.py
                mpirun -np 4 python ex39p.py -ess "Southern Boundary"
                mpirun -np 4 python ex39p.py -src Base

   Description:  This example code demonstrates the use of named attribute
                 sets in MFEM to specify material regions, boundary regions,
                 or source regions by name rather than attribute numbers. It
                 also demonstrates how new named attribute sets may be created
                 from arbitrary groupings of attribute numbers and used as a
                 convenient shorthand to refer to those groupings in other
                 portions of the application or through the command line.
  
                 The particular problem being solved here is nearly the same
                 as that in example 1 i.e. a simple finite element
                 discretization of the Laplace problem -Delta u = 1 with
                 homogeneous Dirichlet boundary conditions and, in this case,
                 an inhomogeneous diffusion coefficient. The diffusion
                 coefficient is given a small default value throughout the
                 domain which is increased by two separate amounts in two named
                 regions.
  
                 This example makes use of a specific input mesh, "compass.msh",
                 containing named domain and boundary regions generated by Gmsh
                 and stored in their "msh" format (version 2.2). This file
                 defines eight boundary regions corresponding to eight compass
                 headings; "ENE", "NNE", "NNW", "WSW", "SSW", "SSE", and "ESE".
                 It also defines nine domain regions; "Base", "N Even", "N Odd",
                 "W Even", "W Odd", "S Even", "S Odd", "E Even", and "E Odd".
                 These regions split the four compass pointers into two halves
                 each and also label the remaining elements as "Base". Starting
                 with these named regions we test the construction of named
                 sets as well as reading and writing these named groupings from
                 and to mesh files.
  
                 The example highlights the use of named attribute sets for
                 both subdomains and boundaries in different contexts as well
                 as basic methods to create named sets from existing attributes.

'''
import os
from os.path import expanduser, join
import numpy as np

import mfem.par as mfem

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)

def run(order=1,
        meshfile='',
        source_name='Rose Even',
        ess_name='Boundary',
        visualization=True):

    # 3. Read the serial mesh from the given mesh file.
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()

    # 4. Refine the serial mesh on all processors to increase the resolution. In
    #    this example we do 'ref_levels' of uniform refinement. We choose
    #    'ref_levels' to be the largest number that gives a final mesh with no
    #    more than 10,000 elements.

    ref_levels = int(np.log(10000./mesh.GetNE())/np.log(2.)/dim)
    for i in range(ref_levels):
        mesh.UniformRefinement()

    # 5. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution. Once the
    #    parallel mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    del mesh

    par_ref_levels = 2
    for l in range(par_ref_levels):
        pmesh.UniformRefinement()

    # 6a. Display attribute set names contained in the initial mesh
    #     GetAttributeSetNames returns Python set object.

    attr_sets = pmesh.attribute_sets
    bdr_attr_sets = pmesh.bdr_attribute_sets

    if myid == 0:
        names = attr_sets.GetAttributeSetNames()
        print("Element Attribute Set Names: " +
              ' '.join(['"' + x + '"' for x in names]))
        names = bdr_attr_sets.GetAttributeSetNames()
        print("Boundary Attribute Set Names: " +
              ' '.join(['"' + x + '"' for x in names]))

    # 6b. Define new regions based on existing attribute sets
    Na = attr_sets.GetAttributeSet("N Even")
    Nb = attr_sets.GetAttributeSet("N Odd")
    Sa = attr_sets.GetAttributeSet("S Even")
    Sb = attr_sets.GetAttributeSet("S Odd")
    Ea = attr_sets.GetAttributeSet("E Even")
    Eb = attr_sets.GetAttributeSet("E Odd")
    Wa = attr_sets.GetAttributeSet("W Even")
    Wb = attr_sets.GetAttributeSet("W Odd")

    # Create a new set spanning the North point
    attr_sets.SetAttributeSet("North", Na)
    attr_sets.AddToAttributeSet("North", Nb)

    # Create a new set spanning the South point
    attr_sets.SetAttributeSet("South", Sa)
    attr_sets.AddToAttributeSet("South", Sb)

    # Create a new set spanning the East point
    attr_sets.SetAttributeSet("East", Ea)
    attr_sets.AddToAttributeSet("East", Eb)

    # Create a new set spanning the West point
    attr_sets.SetAttributeSet("West", Wa)
    attr_sets.AddToAttributeSet("West", Wb)

    # Create a new set consisting of the "a" sides of the compass rose
    attr_sets.SetAttributeSet("Rose Even", Na)
    attr_sets.AddToAttributeSet("Rose Even", Sa)
    attr_sets.AddToAttributeSet("Rose Even", Ea)
    attr_sets.AddToAttributeSet("Rose Even", Wa)

    # Create a new set consisting of the "b" sides of the compass rose
    attr_sets.SetAttributeSet("Rose Odd", Nb)
    attr_sets.AddToAttributeSet("Rose Odd", Sb)
    attr_sets.AddToAttributeSet("Rose Odd", Eb)
    attr_sets.AddToAttributeSet("Rose Odd", Wb)

    # Create a new set consisting of the full compass rose
    Ra = attr_sets.GetAttributeSet("Rose Even")
    Rb = attr_sets.GetAttributeSet("Rose Odd")
    attr_sets.SetAttributeSet("Rose", Ra)
    attr_sets.AddToAttributeSet("Rose", Rb)

    # 6c. Define new boundary regions based on existing boundary attribute sets
    NNE = bdr_attr_sets.GetAttributeSet("NNE")
    NNW = bdr_attr_sets.GetAttributeSet("NNW")
    ENE = bdr_attr_sets.GetAttributeSet("ENE")
    ESE = bdr_attr_sets.GetAttributeSet("ESE")
    SSE = bdr_attr_sets.GetAttributeSet("SSE")
    SSW = bdr_attr_sets.GetAttributeSet("SSW")
    WNW = bdr_attr_sets.GetAttributeSet("WNW")
    WSW = bdr_attr_sets.GetAttributeSet("WSW")

    bdr_attr_sets.SetAttributeSet("Northern Boundary", NNE)
    bdr_attr_sets.AddToAttributeSet("Northern Boundary", NNW)

    bdr_attr_sets.SetAttributeSet("Southern Boundary", SSE)
    bdr_attr_sets.AddToAttributeSet("Southern Boundary", SSW)

    bdr_attr_sets.SetAttributeSet("Eastern Boundary", ENE)
    bdr_attr_sets.AddToAttributeSet("Eastern Boundary", ESE)

    bdr_attr_sets.SetAttributeSet("Western Boundary", WNW)
    bdr_attr_sets.AddToAttributeSet("Western Boundary", WSW)

    bdr_attr_sets.SetAttributeSet("Boundary",
                                  bdr_attr_sets.GetAttributeSet
                                  ("Northern Boundary"))
    bdr_attr_sets.AddToAttributeSet("Boundary",
                                    bdr_attr_sets.GetAttributeSet
                                    ("Southern Boundary"))
    bdr_attr_sets.AddToAttributeSet("Boundary",
                                    bdr_attr_sets.GetAttributeSet
                                    ("Eastern Boundary"))
    bdr_attr_sets.AddToAttributeSet("Boundary",
                                    bdr_attr_sets.GetAttributeSet
                                    ("Western Boundary"))

    # 7. Define a parallel finite element space on the parallel mesh. Here we
    #    use continuous Lagrange finite elements of the specified order. If
    #    order < 1, we instead use an isoparametric/isogeometric space.
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.ParFiniteElementSpace(pmesh, fec)
    size = fespace.GlobalTrueVSize()
    if myid == 0:
        print('Number of finite element unknowns: ' + str(size))

    # 8. Determine the list of true (i.e. parallel conforming) essential
    #    boundary dofs. In this example, the boundary conditions are defined
    #    by marking all the boundary regions corresponding to the boundary
    #    attributes contained in the set named "ess_name" as essential
    #    (Dirichlet) and converting them to a list of true dofs.
    ess_tdof_list = mfem.intArray()
    if bdr_attr_sets.AttributeSetExists(ess_name):
        ess_bdr_marker = bdr_attr_sets.GetAttributeSetMarker(ess_name)
        fespace.GetEssentialTrueDofs(ess_bdr_marker, ess_tdof_list)

    # 9. Set up the parallel linear form b(.) which corresponds to the
    #    right-hand side of the FEM linear system, which in this case is
    #    (1_s,phi_i) where phi_i are the basis functions in fespace and 1_s
    #    is an indicator function equal to 1 on the region defined by the
    #    named set "source_name" and zero elsewhere.

    source_marker = attr_sets.GetAttributeSetMarker(source_name)
    b = mfem.ParLinearForm(fespace)
    one = mfem.ConstantCoefficient(1.0)
    b.AddDomainIntegrator(mfem.DomainLFIntegrator(one), source_marker)
    b.Assemble()

    # 10. Define the solution vector x as a parallel finite element grid
    #     function corresponding to fespace. Initialize x with initial guess of
    #     zero, which satisfies the boundary conditions.
    x = mfem.ParGridFunction(fespace)
    x.Assign(0.0)

    # 11. Set up the parallel bilinear form a(.,.) on the finite element space
    #     corresponding to the Laplacian operator -Delta, by adding the
    #     Diffusion domain integrator.
    a = mfem.ParBilinearForm(fespace)

    defaultCoef = mfem.ConstantCoefficient(1.0e-6)
    baseCoef = mfem.ConstantCoefficient(1.0)
    roseCoef = mfem.ConstantCoefficient(2.0)

    base_marker = attr_sets.GetAttributeSetMarker("Base")
    rose_marker = attr_sets.GetAttributeSetMarker("Rose Even")

    # Impose a very small diffusion coefficient across the entire mesh
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(defaultCoef))

    # Impose an additional, stronger diffusion coefficient in select regions
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(baseCoef), base_marker)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(roseCoef), rose_marker)

    # 12. Assemble the bilinear form and the corresponding linear system,
    #     applying any necessary transformations.
    a.Assemble()

    A = mfem.HypreParMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

    # 13. Solve the system using PCG with hypre's BoomerAMG preconditioner.
    M = mfem.HypreBoomerAMG(A)
    cg = mfem.CGSolver(MPI.COMM_WORLD)
    cg.SetRelTol(1e-12)
    cg.SetMaxIter(2000)
    cg.SetPrintLevel(1)
    cg.SetPreconditioner(M)
    cg.SetOperator(A)
    cg.Mult(B, X)

    # 14. Recover the parallel grid function corresponding to X. This is the
    #     local finite element solution on each processor.
    a.RecoverFEMSolution(X, b, x)

    # 15. Save the refined mesh and the solution in parallel.
    #     This output can be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    #
    #     Python note: In Python, suffix based on myid needs to be added explicitly.
    pmesh.Print('mesh'+smyid)
    x.Save('sol'+smyid)

    # 16. Send the solution by socket to a GLVis server.
    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock << "parallel " << num_procs << " " << myid << "\n"
        sol_sock.precision(8)
        sol_sock << "solution\n" << pmesh << x << "keys Rjmm"
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex1 (Laplace Problem)')
    parser.add_argument('-m', '--mesh',
                        default='compass.msh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help='Finite element order (polynomial degree) or -1 for isoparametric space.')
    parser.add_argument('-src', '--source-attr-name',
                        action='store', default='Rose Even', type=str,
                        help='Name of attribute set containing source.')
    parser.add_argument('-ess', '--ess-attr-name',
                        action='store', default='Boundary', type=str,
                        help='Name of attribute set containing essential BC.')
    parser.add_argument('-no-vis', '--no-visualization',
                        action='store_true',
                        default=False,
                        help='Disable or disable GLVis visualization')

    args = parser.parse_args()
    if myid == 0:
        parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))

    visualization = not args.no_visualization
    source_name = args.source_attr_name
    ess_name = args.ess_attr_name

    run(order=order,
        meshfile=meshfile,
        source_name=source_name,
        ess_name=ess_name,
        visualization=visualization)
