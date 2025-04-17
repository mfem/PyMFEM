'''
   bessel functions  (real order, complex argument)
      iv
      kv
      jv
      hankel1
      hankel2
      yv
'''
import numpy as np
from numpy import exp, log, pi, sin, cos, sqrt, abs, cosh, sinh, sum
from math import gamma
from numba import njit, float64, int64, complex128

# machine precision for double computation
E = np.finfo(np.complex128).resolution  # 1e-15
D = -np.log10(np.finfo(np.complex128).resolution)  # 15

'''
I(n, x)
  n: real
  x: compplex
'''


def iv_smallz(n, z):
    '''
    |z| < 2*sqrt(n+1)
    '''
    a = 1
    sum1 = 0
    for k in range(20):
        sum1 = sum1 + a
        fac = 1/(k+1)/(k+n+1)*(z**2/4)
        a = a*fac
        if abs(fac) < E:
            break

    return sum1*(z/2)**n/gamma(1+n)


def iv_largez(n, z):
    '''
    good for |z| > max(1.2*D+2.4, n**2/2)
    '''
    a = -(4*n**2 - 1)/(8*z)
    b = (4*n**2 - 1)/(8*z)

    sum1 = 1
    sum2 = 1
    for i in range(15):
        k = i+1
        sum1 = sum1 + a
        sum2 = sum2 + b
        k = k+1
        a = -a*(4*n**2 - (2*k-1)**2)/k/(8*z)
        b = b*(4*n**2 - (2*k-1)**2)/k/(8*z)

    ans = exp(z)/sqrt(2*pi*z) * sum1
    if np.imag(z) >= 0:
        ans = ans + sum2*exp(1j*(n+0.5)*pi)*exp(-z)/sqrt(2*pi*z)
    else:
        ans = ans + sum2*exp(-1j*(n+0.5)*pi)*exp(-z)/sqrt(2*pi*z)

    return ans


def iv_middlez(n, z):
    nu = n - np.floor(n)
    size = 40
    data = np.zeros(size, dtype=np.complex128)
    data[-2] = 1

    for i in range(len(data)-2):
        nn = len(data)-2-i+nu
        data[-3-i] = 2*nn/z*data[-2-i] + data[-1-i]

    norm = data[0]
    for k in range(len(data)):
        if k == 0:
            norm = data[0]
            fac = 1.
        else:
            fac = fac*k
            lam = 2*(nu+k)*gamma(k+2*nu)/gamma(1+2*nu)/fac
            norm += data[k]*lam

    data /= norm
    data *= (z/2.)**nu/gamma(1+nu)*exp(z)

    return data[int(np.floor(n))]


def kv_smallz(n, z):
    '''
    |z| > 2, np.real(z) > 0
    '''
    tmp = n-np.floor(n)
    nu = tmp if tmp < 0.5 else tmp-1.

    if nu == -0.5:
        k1 = exp(-z)*sqrt(pi/2/z)
        k2 = exp(-z)*sqrt(pi/2/z)

    if nu != 0:
        g1 = (1./gamma(1-nu) - 1./gamma(1+nu))/(2*nu)
    else:
        nu1 = 1e-5
        g1 = (1/gamma(1-nu1) - 1/gamma(1+nu1))/(2.*nu1)

    g2 = (1/gamma(1-nu) + 1/gamma(1+nu))/2.0

    aa = gamma(1+nu)*gamma(1-nu)

    if abs(nu) == 0.0:
        term1 = 1.0
        term2 = 1.0
    else:
        mu = nu*log(2/z)
        term1 = cosh(mu)
        term2 = sinh(mu)/mu

    fk = aa*(g1*term1 + g2*log(2./z)*term2)
    pk = 0.5*(z/2.)**(-nu)*gamma(1+nu)
    qk = 0.5*(z/2.)**(nu)*gamma(1-nu)
    ck = 1.

    k1 = ck*fk
    k2 = ck*pk
    for i in range(30):
        k = i+1

        fk = (k*fk + pk + qk)/(k**2-nu**2)

        pk = pk/(k-nu)
        qk = qk/(k+nu)
        ck = ck*(z**2/4)/k

        k1 = k1 + ck*fk
        k2 = k2 + ck*(pk - k*fk)

    k2 = 2/z*k2

    num_rec = int(n - nu)

    if num_rec == 0:
        return k1
    if num_rec == 1:
        return k2

    for i in range(num_rec-1):
        nn = i+1+nu
        tmp = k2
        k2 = 2*nn/z*k2+k1
        k1 = tmp

    return k2


def kv_largez(n, z):
    '''
    |z| > 2, np.real(z) > 0
    '''
    tmp = n-np.floor(n)
    nu = tmp if tmp < 0.5 else tmp-1.

    if nu == -0.5:
        k1 = exp(-z)*sqrt(pi/2/z)
        k2 = k1
    else:
        size = 50
        data = np.zeros(size, np.complex128)
        data[1] = 1.

        for i in range(size-2):
            nn = size-i-2.

            an = ((nn-0.5)**2 - nu**2)/nn/(nn+1)
            bn = 2*(nn+z)/(nn+1)

            data[i+2] = -(data[i] - bn*data[i+1])/an

        data /= sum(data)
        data *= (2*z)**(-nu-0.5)

        fac = cos(pi*nu)/pi*gamma(nu+0.5)*gamma(0.5-nu)
        uu0 = data[-1]/fac

        k1 = sqrt(pi)*(2*z)**nu * exp(-z)*uu0
        k2 = k1*(nu+0.5+z - data[-2]/data[-1])/z

    num_rec = int(n - nu)

    if num_rec == 0:
        return k1
    if num_rec == 1:
        return k2

    for i in range(num_rec-1):
        nn = i+1+nu
        tmp = k2
        k2 = 2*nn/z*k2+k1
        k1 = tmp

    return k2


def kv_righthalf(n, z):
    if abs(z) < 2:
        kval = kv_smallz(abs(n), z)
    else:
        kval = kv_largez(abs(n), z)
    return kval


def iv_righthalf(n, z):
    safe_guard = 1.0  # must be >1, expand intermidate region
    if abs(z) < 2*sqrt(abs(n)+1)/safe_guard:
        ival = iv_smallz(abs(n), z)
    elif abs(z) > (1.2*D+2.4)*safe_guard and abs(z) > (abs(n)**2/2)*safe_guard:
        ival = iv_largez(abs(n), z)
    else:
        ival = iv_middlez(abs(n), z)

    if n < 0:
        kval = kv_righthalf(n, z)
        ival = ival + 2./pi*sin(-n*pi)*kval

    if np.imag(z) == 0:
        return np.real(ival)

    return ival


def jv_righthalf(n, z):
    jval = exp(1j*n*pi/2.)*iv_righthalf(n, z/1j)
    if np.imag(z) == 0:
        return np.real(jval)
    return jval

#
#   kv, iv, jv, hankel1, hankel2, yv
#


def kv(n, z):
    if np.real(z) < 0:
        zz = -z
        ival = iv_righthalf(n, zz)
    else:
        zz = z
        ival = 0j

    kval = kv_righthalf(n, zz)

    ee = exp(1j*pi*n)
    if (n-np.floor(n)) == 0.0:
        ee = np.real(ee)

    if np.real(z) < 0 and np.imag(z) >= 0:
        return 1/ee*kval - 1j*pi*ival
    elif np.real(z) < 0 and np.imag(z) < 0:
        return ee*kval + 1j*pi*ival
    else:
        return kval


def iv(n, z):
    if np.real(z) < 0:
        zz = -z
    else:
        zz = z
    ival = iv_righthalf(n, zz)

    ee = exp(1j*pi*n)
    if (n-np.floor(n)) == 0.0:
        ee = np.real(ee)

    if -np.real(z) > 0 and np.imag(z) >= 0:
        return ee*ival
    elif -np.real(z) > 0 and np.imag(z) < 0:
        return 1/ee*ival

    return ival


def jv(n, z):
    if -np.real(z) > 0:
        zz = -z
    else:
        zz = z
    jval = jv_righthalf(n, zz)

    ee = exp(1j*pi*n)
    if (n-np.floor(n)) == 0.0:
        ee = np.real(ee)

    if -np.real(z) > 0 and np.imag(z) >= 0:
        return ee*jval
    elif -np.real(z) > 0 and np.imag(z) < 0:
        return 1/ee*jval

    return jval


def hankel1(n, z):
    zz = z*exp(-1j*pi/2)
    return 2/pi/1j*exp(-1j*pi*n/2)*kv(n, zz)


def hankel2(n, z):
    zz = z*exp(1j*pi/2)
    return -2/pi/1j*exp(1j*pi*n/2)*kv(n, zz)


def yv(n, z):
    return (hankel1(n, z) - hankel2(n, z))/2/1j


jitter = njit(complex128(float64, complex128))

iv_smallz = jitter(iv_smallz)
iv_largez = jitter(iv_largez)
iv_middlez = jitter(iv_middlez)
kv_smallz = jitter(kv_smallz)
kv_largez = jitter(kv_largez)
kv_righthalf = jitter(kv_righthalf)
iv_righthalf = jitter(iv_righthalf)
jv_righthalf = jitter(jv_righthalf)
kv = jitter(kv)
iv = jitter(iv)
jv = jitter(jv)
hankel1 = jitter(hankel1)
hankel2 = jitter(hankel2)
yv = jitter(yv)
