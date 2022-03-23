'''

    Solution of the classical L-shaped domain problem (Poisson)

    William F. Mitchell
    "A collection of 2D elliptic problems for testing adaptive grid refinement algorithms 
    Applied Mathematics and Computation 2013

'''

from math import atan, sqrt
import mfem.ser as mfem
from math import atan2, sqrt, sin, cos
from numba import cfunc, carray
import numpy as np

@cfunc(mfem.scalar_sig)
def LShapedExact(pt, sdim):
    x = pt[0]
    y = pt[1]
    r = sqrt(x*x + y*y)
    alpha = 2. / 3.
    theta = atan2(y, x)
    if x > 0 and abs(y) < 1e-6:
        theta = 0.0
    elif y < 0:
        theta += 2*np.pi
    return r**alpha * sin(alpha * theta)

@cfunc(mfem.vector_sig)
def LShapedExactGrad(pt, out, sdim, vdim):
    out_array = carray(out, (vdim, ))
    x = pt[0]
    y = pt[1]
    alpha = 2. / 3.
    r = sqrt(x*x + y*y)
    if (r == 0):
        r+=1e-12
    theta = atan2(y, x)
    if x > 0 and abs(y) < 1e-6:
        theta = 0.0
    elif y < 0:
        theta += 2*np.pi
    rx = x/r
    ry = y/r
    thetax = - y / r**2
    thetay =   x / r**2
    out_array[0] = alpha * r**(alpha - 1.) *(rx*sin(alpha*theta) + r*thetax * cos(alpha*theta))
    out_array[1] = alpha * r**(alpha - 1.) *(ry*sin(alpha*theta) + r*thetay * cos(alpha*theta))