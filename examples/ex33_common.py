'''
  Ex33_common.py

  This is a translation of MFEM ex33.hpp. LAPACK call, mfem.Vector, 
  mfem.DenseMatrix are replaced using numpy
 

                (Implementation of the AAA algorithm)

  Here, we implement the triple-A algorithm [1] for the rational approximation
  of complex-valued functions,

          p(z)/q(z) ≈ f(z).

  In this file, we always assume f(z) = z^{-α}. The triple-A algorithm
  provides a robust, accurate approximation in rational barycentric form.
  This representation must be transformed into a partial fraction
  representation in order to be used to solve a spectral FPDE.

  More specifically, we first expand the numerator in terms of the zeros of
  the rational approximation,

          p(z) ∝ Π_i (z - z_i),

  and expand the denominator in terms of the poles of the rational
  approximation,

          q(z) ∝ Π_i (z - p_i).

  We then use these zeros and poles to derive the partial fraction expansion

          f(z) ≈ p(z)/q(z) = Σ_i c_i / (z - p_i).

  [1] Nakatsukasa, Y., Sète, O., & Trefethen, L. N. (2018). The AAA algorithm
      for rational approximation. SIAM Journal on Scientific Computing, 40(3),
      A1494-A1522.
'''
import numpy as np
import scipy
from scipy.linalg import eig

import mfem
if mfem.mfem_mode == 'parallel':
    import mfem.par as mfem
    from mfem.par import intArray, doubleArray
else:
    import mfem.ser as mfem
    from mfem.ser import intArray, doubleArray

from sys import float_info
eps = float_info.min


def RationalApproximation_AAA(val, pt, tol, max_order):
    '''
    RationalApproximation_AAA: compute the rational approximation (RA) of data
    val at the set of points pt

    in:
        val        Vector of data values
        pt         Vector of sample points
        tol        Relative tolerance
        max_order  Maximum number of terms (order) of the RA
    out:
        z          Support points of the RA in rational barycentric form
        f          Data values at support points at z
        w          Weights of the RA in rational barycentric form

    See pg. A1501 of Nakatsukasa et al. [1].
    '''

    # number of sample points
    size = len(val)  # .Size()
    assert len(pt) == size, "size mismatch"

    # Initializations
    J = list(range(size))

    c_i = []
    # mean of the value vector
    mean_val = np.mean(val)
    R = np.array([mean_val]*size)

    z = []
    f = []
    w = []

    for k in range(max_order):
        # select next support point
        idx = 0
        tmp_max = 0

        idx = np.argmax(np.abs(val-R))
        # Append support points and data values
        z.append(pt[idx])
        f.append(val[idx])

        # Update index vector
        J.remove(idx)

        # next column in Cauchy matrix
        C_tmp = [(1.0/(pp-pt[idx]) if pp != pt[idx] else np.inf) for pp in pt]

        c_i = np.hstack((c_i, C_tmp))

        h_C = len(C_tmp)
        w_C = k+1

        # note: tranpose is necessary due to the difference of column-major
        # and raw-major of matrix
        C = c_i.reshape(w_C, h_C).transpose()
        Ctemp = C.copy()
        Ctemp = Ctemp*(np.atleast_2d(1/val)).transpose()  # InvLeftScaling
        Ctemp = Ctemp*f                                   # RgithScaling
        A = C - Ctemp
        A = A*(np.atleast_2d(val)).transpose()            # LeftScaling

        h_Am = len(J)
        w_Am = A.shape[1]
        Am = np.zeros((h_Am, w_Am))

        for i in range(h_Am):
            ii = J[i]
            for j in range(w_Am):
                Am[i, j] = A[ii, j]

        u, s, vh = np.linalg.svd(Am)

        w = vh[k, :]

        N = C.dot(w*np.array(f))
        D = C.dot(w)
        R = val.copy()
        for i, ii in enumerate(J):
            R[ii] = N[ii]/D[ii]
        verr = val - R

        if np.max(verr) <= tol*max(val):
            break

    return z, f, w


def ComputePolesAndZeros(z, f, w):
    '''
    ComputePolesAndZeros: compute the poles  and zeros of the
    rational function f(z) = C p(z)/q(z) from its ration barycentric form.

    in:
        z      Support points in rational barycentric form
        f      Data values at support points @a z
        w      Weights in rational barycentric form
    out:
        poles  Array of poles (roots of p(z))
        zeros  Array of zeros (roots of q(z))
        scale  Scaling constant in f(z) = C p(z)/q(z)

    See pg. A1501 of Nakatsukasa et al. [1].
    '''

    # Initialization
    poles = []
    zeros = []

    # Compute the poles
    m = len(w)
    B = np.zeros((m+1, m+1))
    E = np.zeros((m+1, m+1))

    for i in range(m+1):
        if i == 0:
            continue
        B[i, i] = 1.
        E[0, i] = w[i-1]
        E[i, 0] = 1.
        E[i, i] = z[i-1]

    # real part of eigen value
    evalues = eig(E, B, left=False, right=False).real
    new_poles = evalues[np.isfinite(evalues)]

    poles.extend(new_poles)

    B = np.zeros((m+1, m+1))
    E = np.zeros((m+1, m+1))
    for i in range(m+1):
        if i == 0:
            continue
        B[i, i] = 1.
        E[0, i] = w[i-1] * f[i-1]
        E[i, 0] = 1.
        E[i, i] = z[i-1]

    # real part of eigen value
    evalues = eig(E, B, left=False, right=False).real
    new_zeros = evalues[np.isfinite(evalues)]

    zeros.extend(new_zeros)

    scale = np.dot(w, f)/np.sum(w)

    return poles, zeros, scale


def PartialFractionExpansion(scale, poles, zeros):
    '''
    PartialFractionExpansion: compute the partial fraction expansion of the
    rational function f(z) = Σ_i c_i / (z - p_i) from its poles and zeros
    @a zeros [in].

    in: 
        poles   Array of poles (same as p_i above)
        zeros   Array of zeros
        scale   Scaling constant
    out:
        coeffs  Coefficients c_i 
    '''

    # Note: C p(z)/q(z) = Σ_i c_i / (z - p_i) results in an system of equations
    # where the N unknowns are the coefficients c_i. After multiplying the
    # system with q(z), the coefficients c_i can be computed analytically by
    # choosing N values for z. Choosing z_j = = p_j diagonalizes the system and
    # one can obtain an analytic form for the c_i coefficients. The result is
    # implemented in the code block below.

    psize = len(poles)
    zsize = len(zeros)
    coeffs = [scale] * psize

    for i in range(psize):
        tmp_numer = 1.0
        for j in range(zsize):
            tmp_numer *= poles[i]-zeros[j]

        tmp_denom = 1.0
        for k in range(psize):
            if k != i:
                tmp_denom *= poles[i]-poles[k]

        coeffs[i] *= tmp_numer / tmp_denom
    return coeffs


def ComputePartialFractionApproximation(alpha,
                                        lmax=1000.,
                                        tol=1e-10,
                                        npoints=1000,
                                        max_order=100):
    '''
    ComputePartialFractionApproximation: compute a rational approximation (RA)
    in partial fraction form, e.g., f(z) ≈ Σ_i c_i / (z - p_i), from sampled
    values of the function f(z) = z^{-a}, 0 < a < 1.

    in:
       alpha         Exponent a in f(z) = z^-a
       lmax,npoints  f(z) is uniformly sampled @a npoints times in the
                     interval [ 0, @a lmax ]
       tol           Relative tolerance
       max_order     Maximum number of terms (order) of the RA

    out:
       coeffs        Coefficients c_i
       poles         Poles p_i
    '''

    assert alpha < 1., "alpha must be less than 1"
    assert alpha > 0., "alpha must be greater than 0"
    assert npoints > 2, "npoints must be greater than 2"
    assert lmax > 0,  "lmin must be greater than 0"
    assert tol > 0,  "tol must be greater than 0"

    dx = lmax / (npoints-1)

    x = np.arange(npoints)*dx
    val = x**(1-alpha)

    # Apply triple-A algorithm to f(x) = x^{1-a}
    z, f, w = RationalApproximation_AAA(val,  # mfem.Vector(val),
                                        x,  # mfem.Vector(x),
                                        tol, max_order)

    # Compute poles and zeros for RA of f(x) = x^{1-a}
    poles, zeros, scale = ComputePolesAndZeros(z, f, w)

    # Remove the zero at x=0, thus, delivering a RA for f(x) = x^{-a}
    zeros.remove(0.0)

    # Compute partial fraction approximation of f(x) = x^{-a}
    coeffs = PartialFractionExpansion(scale, poles, zeros)

    return poles, coeffs
