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

def run(mesh_file="",
        order=1,
        ser_ref_levels=0,
        par_ref_levels=0,
        visualization=1,
        prob=0,
        mode=1,
        herm_conv=True,
        a_coef=0.0,
        epsilon=1.0,
        sigma=20.,
        mu=1.0,
        device='cpu',
        pa=False,
        port_bc_attr = None,
        freq=-1.0,):

    mesh_file = "../data/fichera-mixed.mesh"
    if not mixed or pa:
        mesh_file = "../data/fichera.mesh"
        
    if a_coef != 0.0:
         mu_ = 1.0 / a_coef
    omega = 2*pi
    if freq > 0.0:
        omega = 2.0 * pi * freq

    if (len(port_bc_attr) == 0 and
        (mesh_file == "../data/fichera-mixed.mesh" or
         mesh_file == "../data/fichera.mesh")):
        port_bc_attr = mfem.intArray([7, 8, 11, 12])

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
        pmesh.UniformRefinement()

    # 6a. Define a parallel mesh by a partitioning of the serial mesh. Refine
    #    this mesh further in parallel to increase the resolution. Once the
    #    parallel mesh is defined, the serial mesh can be deleted.
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    for i in range(par_ref_levels):
        pmesh.UniformRefinement()
    
    # 6b. Extract a submesh covering a portion of the boundary
    pmesh_cond = mfem.ParrSubMesh.CreateFromBoundary(pmesh, port_bc_attr)

        
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
                       help="Damping coefficient (or sigma).");
    parser.add_argument("-mu", "--permeability",
                        action='store', type=float, default=1.0,
                        help="Permeability of free space (or 1/(spring constant)).")
    parser.add_argument("-eps", "--permittivity",
                        action='store', type=float, default=1.0,
                        help="Permittivity of free space (or mass constant).")
    parser.add_argument("-sigma", "--conductivity",
                        action='store', type=float, default=20.0,
                        help="Conductivity (or damping constant).")
    parser.add_argument("-f", "--frequency",
                        action='store',
                        type=float,
                        default=-1.0,
                        help="Set the frequency for the exact")
    parser.add_argument("-pbc", "--port-bc-attr",
                        action='store', type=str, default="",
                       "Attributes of port boundary condition")
    parser.add_argument("-herm", "--hermitian",
                        action='store_true',
                        default=True,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument("-no-herm", "--no-hermitian",
                        action='store_true',
                        default=False,
                        help="Do not use convention for Hermitian operators.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=True,
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

    port_bc_attr = [int(x) for x in args.port_bc_attr.split(' ')]
    globals()["mixed"] = not args.hex_mesh
    
    run(mesh_file=meshfile,
        order=args.order,
        ser_ref_levels=args.refine_serial,
        par_ref_levels=args.refine_parallel,
        prob=args.problem_type,
        mode=args.eigenmode,
        visualization=args.visualization,
        # visualization=False,
        herm_conv=herm,
        a_coef=args.stiffness_coef,
        epsilon=args.permittivity,
        sigma=args.conductivity,
        mu=args.permeability,
        device=args.device,
        pa=args.partial_assembly,
        port_bc_attr=port_bc_attr,
        freq=args.frequency)

