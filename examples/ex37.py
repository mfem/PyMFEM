'''
   PyMFEM example 37

   See c++ version in the MFEM library for more detail

   Sample runs:
       python ex37.py -alpha 10
       python ex37.py -alpha 10 -pv
       python ex37.py -lambda 0.1 -mu 0.1
       python ex37.py -o 2 -alpha 5.0 -mi 50 -vf 0.4 -ntol 1e-5
       python ex37.py -r 6 -o 1 -alpha 25.0 -epsilon 0.02 -mi 50 -ntol 1e-5


'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor

visualization = True
paraview_output = False


def run(ref_levels=5,
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
    pass


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
    parser.add_argument("-lambda", "--lambda",
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

    run(ref_levels=args.refine,
        order=args.order,
        alpha=args.alpha,
        epsilon=args.epsilon_thickness,
        vol_fraction=args.volume_fraction,
        max_it=args.max_it,
        itol=args.abs_tol,
        ntol=args.ref_tol,
        rho_min=args.psi_min,
        llambda=args.lambda,
        mu=args.mu)
        use_umfpack = use_umfpack)
