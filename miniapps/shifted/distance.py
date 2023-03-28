'''
   Distance.py
      See c++ version in the MFEM library for more detail 
'''
import os
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser, join, dirname
import numpy as np
from numpy import sin, cos, exp, sqrt, pi, abs, array, floor, log

from mpi4py import MPI
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '.'+'{:0>6d}'.format(myid)


def VisualizeField(vishost,
                   visport,
                   gf,
                   title,
                   x, y, w, h, keys=""):
    
     mesh = gf.FESpace().GetMesh()

     newly_opened = False
     connection_failed = 0

     make_connection = False

     sock = mfem.socketstream(vishost, visport)
     sock.precision(8)

     sock << "solution\n"
     sock << mesh << gf
     
     sock << "window_title '" << title << "'\n"
     sock << "window_geometry "
     sock << x << " " << y << " " << w << " " << h << "\n"
     if keys != "":
         sock << "keys " << keys << "\n"
     else:
         sock << "keys maaAc\n"
     sock.endline()

def run(order=2,
        problem=1,
        solver_type=0,
        rs_levels=2,
        mesh_file="",
        t_param=1.0,
        visualization=True,
        device='cpu',
        numba=True):

    device = mfem.Device(device)
    device.Print()

    mesh = mfem.Mesh(mesh_file, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()    
    for i in range(rs_levels):
        mesh.UniformRefinement()

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    mesh.Clear()

    if problem == 0:
        ls_coeff = mfem.DeltaCoefficient(0.5, -0.5, 1000.0)
        smooth_steps = 0
        
    elif problem == 1:
        if numba:
            @mfem.jit.scalar()
            def ls_coeff(x):
                if len(x) == 2:
                     xc = x[0] - 0.5
                     yc = x[1] - 0.5
                     r = sqrt(xc*xc + yc*yc)
                     return -1.0 if r >= 0.4 else 1.0
                elif len(x) == 3:
                     xc = x[0] - 0.
                     yc = x[1] - 0.
                     zc = x[2] - 0.                 
                     r = sqrt(xc*xc + yc*yc + zc*zc)
                     return -1.0 if r >= 0.8 else 1.0
                else:
                     return -1.0 if x[0] >= 0.5 else 1.0
        else:
            assert False, "distance requires numba"
        smooth_steps = 0
        
    elif problem == 2:
        if numba:
            @mfem.jit.scalar()
            def ls_coeff(x):
                sine = (0.25 * sin(4 * pi * x[0]) +
                        0.05 * sin(16 * pi * x[0]))
                return -1.0 if x[1] >= sine + 0.5 else 1.0
        else:
            assert False, "distance requires numba"
        smooth_steps = 0
        
    elif problem == 3:
        if numba:
            @mfem.jit.scalar()
            def ls_coeff(xx):
                period = 2.0 * pi
                x=xx[0]*period
                y=xx[1]*period
                z = 0.0
                if len(xx) == 3:
                    z=xx[2]*period
                return (sin(x)*cos(y) + sin(y)*cos(z) + sin(z)*cos(x))
        else:
            assert False, "distance requires numba"
        smooth_steps = 0
                        
    dx = mfem.dist_solver.AvgElementSize(pmesh)

    if solver_type == 0:
        ds = mfem.dist_solver.HeatDistanceSolver(t_param * dx * dx)
        if problem == 0:
            ds.transform = False
        ds.mooth_steps = smooth_steps
        ds.vis_glvis = False
        dist_solver = ds
    elif solver_type == 1:
        p = 10
        newton_iter = 50
        ds = mfem.dist_solver.PLapDistanceSolver(p, newton_iter)
        dist_solver = ds
    else:
         assert False, "Wrong solver option."

    fec = mfem.H1_FECollection(order, dim)
    pfes_s = mfem.ParFiniteElementSpace(pmesh, fec)
    pfes_v = mfem.ParFiniteElementSpace(pmesh, fec, dim)
                        
    distance_s = mfem.ParGridFunction(pfes_s)
    distance_v = mfem.ParGridFunction(pfes_v)

    # Smooth-out Gibbs oscillations from the input level set. The smoothing
    # parameter here is specified to be mesh dependent with length scale dx.
    filt_gf = mfem.ParGridFunction(pfes_s)
    filter = mfem.dist_solver.PDEFilter(pmesh, 1.0 * dx)
                        
    if problem != 0:
        filter.Filter(ls_coeff, filt_gf)
    else:
        filt_gf.ProjectCoefficient(ls_coeff)

    ls_filt_coeff = mfem.GridFunctionCoefficient(filt_gf)

    dist_solver.ComputeScalarDistance(ls_filt_coeff, distance_s)
    dist_solver.ComputeVectorDistance(ls_filt_coeff, distance_v)

    if visualization:
        size = 500
        vishost = "localhost"
        visport = 19916
                        
        VisualizeField(vishost, visport, filt_gf,
                       "Input Level Set", 0, 0, size, size);

        MPI.COMM_WORLD.Barrier()
                        
        VisualizeField(vishost, visport, distance_s,
                       "Distance", size, 0, size, size,
                       "rRjmm********A");

        MPI.COMM_WORLD.Barrier()
                        
        VisualizeField(vishost, visport, distance_v,
                       "Distances", 2*size, 0, size, size,
                       "rRjmm********vveA");

    zero = mfem.ConstantCoefficient(0.0)
    d_norm  = distance_s.ComputeL2Error(zero)
    if myid == 0:
       print("Norm: " +  "{:g}".format(d_norm))

if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Disntance (translated from miniapps/shifted/distance.cpp)')

    parser.add_argument('-m', '--mesh',
                        default='inline-quad.mesh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-s', '--solver',
                        action='store', type=int, default=0,
                        help= "\n".join(["Solver type",
                                         "\t 0: Heat"
                                         "\t 1: P-Laplacian"]))
    parser.add_argument('-p', '--problem-type',
                        action='store', type=int, default=1,
                        help="\n".join(["Problem type:",
                                        "\t 0: Point source",
                                        "\t 1: Circle / sphere level set in 2D / 3D",
                                        "\t 2: 2D sine-looking level set",
                                        "\t 3: Gyroid level set in 2D or 3D"]))
    parser.add_argument('-o', '--order',
                        action='store', default=2, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-rs', '--refine-serial',
                        action='store', default=2, type=int,
                        help="Number of times to refine the mesh uniformly in serial")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument("-t", "--t-param",
                        default=1., type=float,
                        help="Diffusion time step (scaled internally scaled by dx*dx).")
    parser.add_argument("-n", "--numba",
                        default=1, action='store', type=int,
                        help="Use Number compiled coefficient")
    parser.add_argument("-d", "--device",
                        default="cpu", type=str,
                        help="Device configuration string, see Device::Configure().")
    
    args = parser.parse_args()
    parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', '..', 'data', args.mesh))
    numba = (args.numba == 1)

    run(order=args.order,
        problem=args.problem_type,
        solver_type=args.solver,
        rs_levels=args.refine_serial,
        mesh_file=meshfile,
        #visualization=args.visualization,
        visualization=True,
        t_param=args.t_param,
        numba=numba,
        device=args.device)

    
