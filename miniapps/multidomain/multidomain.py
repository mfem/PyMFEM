'''
   multidomain.py

   See c++ version in the MFEM library for more detail
'''
from mpi4py import MPI
import mfem.par as mfem
from mfem.par import intArray, doubleArray
import os
import sys
from os.path import expanduser, join


num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)

def run(dt=1e-5,
        t_final=5,
        order=2,
        visualization=True,
        vis_step=10):
    

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
    fes_cylinder = mfem.ParFiniteElementSpace (cylinder_submesh, fec)
    
    inflow_attributes = intArray([0]*
        cylinder_submesh.bdr_attributes.Max())
    inflow_attributes[7] = 1

    inner_cylinder_wall_attributes = intArray([0]*
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
                                     interface_tdofs);
    ess_tdofs.Append(inflow_tdofs)
    ess_tdofs.Append(interface_tdofs)
    ess_tdofs.Sort()
    ess_tdofs.Unique()
    
    cd_tdo = ConvectionDiffusionTDO(fes_cylinder, ess_tdofs);

   ParGridFunction temperature_cylinder_gf(&fes_cylinder);
   temperature_cylinder_gf = 0.0;
    '''

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


    try:
        from numba import jit
        HAS_NUMBA = True
    except ImportError:
        assert False, "This example requires numba to run"
    
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
        vis_step=args.visualization_steps)

    
