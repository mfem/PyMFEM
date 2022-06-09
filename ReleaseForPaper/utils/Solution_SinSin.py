'''

    Solution of the classical L-shaped domain problem (Poisson)

    William F. Mitchell
    "A collection of 2D elliptic problems for testing adaptive grid refinement algorithms 
    Applied Mathematics and Computation 2013

'''

import mfem.ser as mfem
from math import sin, cos
from numba import cfunc, carray
import numpy as np

@cfunc(mfem.scalar_sig)
def SinSinExact(pt, sdim):
    x = pt[0]
    y = pt[1]
    return sin(np.pi*x) * sin(np.pi*y)

@cfunc(mfem.vector_sig)
def SinSinExactGrad(pt, out, sdim, vdim):
    out_array = carray(out, (vdim, ))
    x = pt[0]
    y = pt[1]
    out_array[0] = np.pi * cos(np.pi*x) * sin(np.pi*y)
    out_array[1] = np.pi * sin(np.pi*x) * cos(np.pi*y)

@cfunc(mfem.scalar_sig)
def SinSinExactLaplace(pt, sdim):
    x = pt[0]
    y = pt[1]
    return 2 * np.pi * np.pi * sin(np.pi*x) * sin(np.pi*y)