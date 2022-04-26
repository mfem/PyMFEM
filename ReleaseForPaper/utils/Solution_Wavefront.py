'''

    Wavefront parameterization from William F. Mitchell
    "A collection of 2D elliptic problems for testing adaptive grid refinement algorithms 
    Applied Mathematics and Computation 2013

'''

from math import atan, sqrt
import mfem.ser as mfem
import numpy as np

class WavefrontSolutionCoefficient(mfem.PyCoefficient):
    def __init__(self, solntype='steep'):
        mfem.PyCoefficient.__init__(self)
        self.solntype = solntype
        if self.solntype == 'mild':
            alpha = 20.0
            xc = -0.05
            yc = -0.05
            r0 = 0.7
        elif self.solntype == 'steep':
            alpha = 1000.0
            xc = -0.05
            yc = -0.05
            r0 = 0.7
        elif self.solntype == 'asymmetric':
            alpha = 1000.0
            xc = 1.5
            yc = 0.25
            r0 = 0.92
        elif self.solntype == 'well':
            alpha = 50.0
            xc = 0.5
            yc = 0.5
            r0 = 0.25
        else:
            assert(True, 'unknown solution type')
        self.alpha = alpha
        self.xc = xc
        self.yc = yc
        self.r0 = r0

    def EvalValue(self, x):
        alpha = self.alpha
        xc = self.xc
        yc = self.yc
        r0 = self.r0
        r = sqrt((x[0] - xc)**2 + (x[1] - yc)**2)
        return atan(alpha * (r - r0))

class WavefrontRHSCoefficient(WavefrontSolutionCoefficient):
    def __init__(self, solntype='steep'):
        super().__init__(solntype)
    def EvalValue(self, x):
        alpha = self.alpha
        xc = self.xc
        yc = self.yc
        r0 = self.r0
        r = sqrt((x[0] - xc)**2 + (x[1] - yc)**2)
        num = - ( alpha - alpha**3 * (r**2 - r0**2) )
        denom = r * ( alpha**2 * r0**2 + alpha**2 * r**2 - 2 * alpha**2 * r0 * r + 1.0 )**2
        denom = max(denom,1e-8)
        return num / denom