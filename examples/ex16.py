'''
   MFEM example 16

   How to run:
      mpirun -np 2 python2.7 <arguments>

   Example of arguments:
      ex16.py -m inline-tri.mesh
      ex16.py -m disc-nurbs.mesh -tf 2
      ex16.py -s 1 -a 0.0 -k 1.0
      ex16.py -s 2 -a 1.0 -k 0.0
      ex16.py -s 3 -a 0.5 -k 0.5 -o 4
      ex16.py -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
      ex16.py -m fichera-q2.mesh
      ex16.py -m escher.mesh
      ex16.py -m beam-tet.mesh -tf 10 -dt 0.1
      ex16.py -m amr-quad.mesh -o 4 -r 0
      ex16.py -m amr-hex.mesh -o 2 -r 0

'''
import sys
import argparse
from os.path import expanduser, join
import numpy as np
from mfem import path

import mfem.ser as mfem
from mfem.ser import intArray

class ConductionOperator(mfem.PyTimeDependentOperator):
    def __init__
    def Mult(self, u, u_dt):    
    def ImplicitSolve(self, dt, u, k):
    def SetParameters(self, u):

class InitialTemperature(mfem.PyCoefficient):
    def EvalValue(self, x):
        xx = np.array(x)
        norm2 = float(np.sum(xx**2) )
        if norm2 < 0.5: return 2.0
        return 1.0

parser = argparse.ArgumentParser(description='Ex17')
parser.add_argument('-m', '--mesh',
                    default = 'star.mesh',
                    action = 'store', type = str,
                    help='Mesh file to use.')
parser.add_argument('-r', '--refine',
                    action = 'store', default = -1, type=int,
       help = "Number of times to refine the mesh uniformly, -1 for auto.")
parser.add_argument('-o', '--order',
                    action = 'store', default = 1, type=int,
                    help = "Finite element order (polynomial degree)");
parser.add_argument('-s', '--ode-solver',
                    action = 'store' default = 3, type = int, 
      help = '\n'.join(["ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3",
                "\t\t 11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4."]))
parser.add_argument('-t', '--t-final',
                    action = 'store', default = 0.5, type=float,
                    help = "Final time; start time is 0.")
parser.add_argument('-a', '--alpha',
                    action = 'store', default = 0.01, type=float,
                    help = 'Alpha coefficient')
parser.add_argument('-k', '--kappa',
                    action = 'store', default = 0.5, type=float,
                    help = 'Kappa coefficient')
parser.add_argument('-vis', '--visualization',
                    action = 'store_true',
                    help='Enable GLVis visualization')
parser.add_argument('-vs', '--visualization-steps',
                    action = 'store', default = 5, type = int,  
                    help = "Visualize every n-th timestep.");                    

args = parser.parse_args()
ref_levels = args.refine
order = args.order
alpha = args.alpha;
kappa = args.kappa;
visualization = args.visualization
vis_step = args.visualization_step

