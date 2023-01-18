'''
   MFEM example 33

   See c++ version in the MFEM library for more detail

   Sample runs:  python ex33.py -m ../data/square-disc.mesh -alpha 0.33 -o 2
                 python ex33.py -m ../data/square-disc.mesh -alpha 4.5 -o 3
                 python ex33.py -m ../data/star.mesh -alpha 1.4 -o 3
                 python ex33.py -m ../data/star.mesh -alpha 0.99 -o 3
                 python ex33.py -m ../data/inline-quad.mesh -alpha 0.5 -o 3
                 python ex33.py -m ../data/amr-quad.mesh -alpha 1.5 -o 3
                 python ex33.py -m ../data/disc-nurbs.mesh -alpha 0.33 -o 3
                 python ex33.py -m ../data/disc-nurbs.mesh -alpha 2.4 -o 3 -r 4
                 python ex33.py -m ../data/l-shape.mesh -alpha 0.33 -o 3 -r 4
                 python ex33.py -m ../data/l-shape.mesh -alpha 1.7 -o 3 -r 5
    Verification runs:
                 python ex33.py -m ../data/inline-segment.mesh -ver -alpha 1.7 -o 2 -r 2
                 python ex33.py -m ../data/inline-quad.mesh -ver -alpha 1.2 -o 2 -r 2
                 python ex33.py -m ../data/amr-quad.mesh -ver -alpha 2.6 -o 2 -r 2
                 python ex33.py -m ../data/inline-hex.mesh -ver -alpha 0.3 -o 2 -r 1
'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor

def run(alpha=0.5,
        order=1,
        rs=3,
        meshfile=None,
        visualization=False,
        vef=False):


   coeffs = doubleArray()
   poles  = doubleArray()
   progress_steps = 1;

   power_of_laplace = floor(alpha)
   exponent_to_approximate = alpha - power_of_laplace

   if exponent_to_approximate > 1e-12:
       print("Approximating the fractional exponent " + 
             str(exponent_to_approximate))
       ComputePartialFractionApproximation(exponent_to_approximate,
                                          coeffs,
                                          poles)

       # If the example is build without LAPACK, the exponent_to_approximate
       # might be modified by the function call above.
       alpha = exponent_to_approximate + power_of_laplace;
       integer_order = False
   else:
       print("Treating integer order PDE.")       
       integer_order = True;

   



    
if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex33 (Fractional PDE)')
    parser.add_argument('-m', '--mesh',
                        default="star.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space")
    parser.add_argument('-r', '--refs',
                        action='store', default=3, type=int,
                        help="Number of times to refine the mesh uniformly in serial")
    parser.add_argument("-alpha", "--alpha",
                        action='store',
                        type=float,
                        default=0.5,
                        help="Fractional exponent")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument('-ver', '--verification',
                        action='store_true',
                        help="Use sinusoidal function (f) for analytic comparison.")

    try:
        from numba import jit
        HAS_NUMBA = True
    except ImportError:
        assert False, "This example requires numba to run"
    
    args = parser.parse_args()
    parser.print_options(args)

    order = args.order
    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    visualization = args.visualization
    alpha = args.alpha


    run(alpha=alpha,
        order=order,
        rs=args.refs,
        meshfile=meshfile,
        visualization=visualization,
        vef=args.verifification)

    

    
