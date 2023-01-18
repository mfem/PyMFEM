'''
  Ex33_common.py

  This is a line-by-line translation of MFEM ex33.hpp. MFEM_USE_LAPACK is skipped

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

import mfem
if mfem.mfem_mode == 'parallel':
    import mfem.par as mfem
else:
    import mfem.ser as mfem

from mfem import intArray, doubleArray    

def RationalApproximation_AAA(val, pt, z, f, w, tol, max_order):
    # number of sample points
    size = val.Size()
    assert pt.Size() == size, "size mismatch"

    # Initializations
    J = intArray(list(range(size)))
    z.SetSize(0)
    f.SetSize(0)


    c_i = doubleArray()
    C = mfem.DenseMatrix()
    Ctemp = mfem.DenseMatrix()    
    A = mfem.DenseMatrix()
    Am = mfem.DenseMatrix()
    f_vec = mfem.Vector()
    
    # mean of the value vector    
    mean_val = val.Sum()/size
    R = mfem.Vector([mean_val]*val.Size())

    for k in range(max_order):
        # select next support point
        idx = 0
        tmp_max = 0
        
        for j in range(size):
            tmp = abs(val[j] - R[j])
            if tmp > tmp_max:
               tmp_max = tmp
               idx = j

        # Append support points and data values
        z.Append(pt[idx])
        f.Append(val[idx])

        # Update index vector
        J.DeleteFirst(idx);

        # next column in Cauchy matrix
        C_tmp = doubleArray(size);
        for j in range(size):
            C_tmp[j] = 1.0/(pt[j]-pt[idx])
            
        c_i.Append(C_tmp);
        h_C = C_tmp.Size()
        w_C = k+1
        
        C.UseExternalData(c_i.GetData(), h_C, w_C)
        Ctemp.Assign(C)

        f_vec.SetDataAndSize(f.GetData(),f.Size());        

        Ctemp.InvLeftScaling(val)
        Ctemp.RightScaling(f_vec)

        A.SetSize(C.Height(), C.Width())
        mfem.Add(C, Ctemp, -1.0, A)
        A.LeftScaling(val)

        h_Am = J.Size()
        w_Am = A.Width()
        Am.SetSize(h_Am,w_Am)
        
        for i in range(h_Am):
            ii = J[i];
            for j in range(w_Am):
                Am[i,j] = A[ii,j]


        AMM = Am.GetDataArray()
        u, s, vh = np.linalg.svd(AMM, full_matrices=True)

        w.Assign(v[k,:])

       # N = C*(w.*f); D = C*w; % numerator and denominator
       aux = mfem.Vector(w)
       aux *= f_vec;
       N = mfem.Vector(C.Height()) # Numerator
       C.Mult(aux, N)
       D =mfem.Vector(C.Height())  # Denominator
       C.Mult(w,D);

       R.Assign(val)
       for i in range(J.Size()):
           ii = J[i];
           R[ii] = N[ii]/D[ii]
       }

       val = mfem.Vector(val)
       verr -= R

       if (verr.Normlinf() <= tol*val.Normlinf()):
           break
   
    
def ComputePolesAndZeros(z, f, w, poles, zeros, scale):
    
    # Initialization
    poles.SetSize(0);
    zeros.SetSize(0);

    # Compute the poles
    m = w.Size();
    B = np.zeros((m+1, m+1))
    E = np.zeros((m+1, m+1))
    
    for i in range(m):
        B[i,i] = 1.
        E[0,i] = w[i-1]
        E[i,0] = 1.
        E[i,i] = z[i-1]

    evalues = scipy.linalg.eig(E, B)
    new_poles = evalues[np.isfinite(evalues)]

    for x in new_poles:
        poles.Append(x)
    

    scale = w*f/w.Sum()
    


