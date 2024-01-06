'''
   PyMFEM example 37 - Serial/Parallel Shared code
'''

import os
from os.path import expanduser, join

from numpy import sqrt, log, exp
import numpy as np
from numba import njit
from numba.types import float64

from mfem import mfem_mode

if mfem_mode == 'serial':
    import mfem.ser as mfem
    from mfem.ser import intArray, doubleArray
    use_parallel = False
else:
    import mfem.par as mfem
    from mfem.par import intArray, doubleArray
    use_parallel = True


@njit(float64(float64))
def inv_sigmoid(x):
    '''
    Inverse sigmoid function
    '''
    tol = 1e-12
    x = min(max(tol, x), 1.0-tol)
    return log(x/(1.0-x))


@njit(float64(float64))
def sigmoid(x):
    '''
    Sigmoid function
    '''
    if x >= 0:
        return 1.0/(1.0 + exp(-x))
    else:
        return exp(x)/(1.0 + exp(x))


@njit(float64(float64))
def der_sigmoid(x):
    '''
    Derivative of sigmoid function    
    '''
    tmp = sigmoid(-x)
    return tmp - tmp**2


def MappedGridFunctionCoefficient(gf, func):

    c_gf = mfem.GridFunctionCoefficient(gf)

    @mfem.jit.scalar(dependency=(c_gf,))
    def coeff(ptx, c_gf):
        return func(c_gf)
    return coeff


def DiffMappedGridFunctionCoefficient(gf, other_gf, func, comp=1):
    c_gf = mfem.GridFunctionCoefficient(gf)
    c_ogf = mfem.GridFunctionCoefficient(other_gf)

    @mfem.jit.scalar(dependency=(c_gf, c_ogf))
    def coeff(ptx, c_gf, c_ogf):
        return func(c_gf) - func(c_ogf)
    return coeff


class SIMPInterpolationCoefficient():
    '''
    Python Note: 
       Approach here is to replace Eval in C++ example using the dependency 
       feature of mfem.jit.

       In order to avoid repeating Numba-Jitting in iteration loop, we use
       SetGridFunction to update the GridFunction referred from 
       GridFunctionCoefficient.
    '''

    def __init__(self, rho_filter, min_val=1e-6, max_val=1.0, exponent=3):
        val = mfem.GridFunctionCoefficient(rho_filter)

        @mfem.jit.scalar(dependency=(val,))
        def coeff(ptx, val):
            coeff = min_val + val**exponent*(max_val-min_val)
            return coeff

        self.c_gf = val
        self.coeff = coeff

    def Update(self, rho_filter):
        self.c_gf.SetGridFunction(rho_filter)


def VolumeForceCoefficient(r, center, force):

    @mfem.jit.vector(shape=(len(center),))
    def coeff(ptx):
        cr = sqrt(sum((ptx - center)**2))
        if cr < r:
            return np.array((force[0], force[1]))
        else:
            return np.array((0.0, 0.0))
    return coeff
