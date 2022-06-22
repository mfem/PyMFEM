'''
    Python class for computing element statistics
'''

from collections import namedtuple
import numpy as np
from math import *
from scipy.stats import describe

stats = namedtuple('Statistics',
                ('nels', 'min', 'max', 'mean', 'variance', 'skewness', 'kurtosis'))

########################
# local error statistics
########################
def Statistics(eta, num_dofs, p=None, d=2):  # prev num_dofs = np.e; p=order; d=dimension 2 or 3
    if p is None:
        zeta = np.log(len(eta)*np.abs(eta)**2)/np.log(num_dofs)
        # zeta = np.log(len(eta)*np.abs(eta)**2)/np.log(len(eta))
    else:
        zeta = np.sqrt(len(eta)) * np.power(num_dofs, p/d) * eta
    description = describe(zeta,bias=False)
    nels = description.nobs
    min = description.minmax[0]
    max = description.minmax[1]
    mean = description.mean
    variance = description.variance
    skewness = description.skewness
    kurtosis = description.kurtosis
    return stats(nels, min, max, mean, variance, skewness, kurtosis)

########################
# total estimated error
########################
def GlobalError(eta):
    return np.sqrt(np.sum(np.abs(eta)**2))