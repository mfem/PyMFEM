'''
   PyMFEM example 38

   See c++ version in the MFEM library for more detail

   Sample runs:
       python ex38.py
       python ex38.py -i volumetric1d
       python ex38.py -i surface2d
       python ex38.py -i surface2d -o 4 -r 5
       python ex38.py -i volumetric2d
       python ex38.py -i volumetric2d -o 4 -r 5
       python ex38.py -i surface3d
       python ex38.py -i surface3d -o 4 -r 5
       python ex38.py -i volumetric3d
       python ex38.py -i volumetric3d -o 4 -r 5

'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray

import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor

visualization = True


def run(ref_levels=3,
        order=2,
        inttype="surface2d"):
    pass 

if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser
    
    parser.add_argument("-o", "--order",
                        action='store', default=2, type=int,
                        help="Order (degree) of the finite elements.")
    parser.add_argument("-r", "--refine",
                        action='store', default=3, type=int,
                        help="Number of meh refinements")
    parser.add_argument(inttype, "-i", "--integrationtype",
                        action='store', default="surface2d", type=string,                        
                        help="IntegrationType to demonstrate");
    parser.add_argument("-no-vis", "--no-visualization",
                        action='store_true', default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    globals()["visualization"] = not args.no_visualization
    
    run(ref_levels=args.refine,
        order=args.order,
        inttype=args.integrationtype)

