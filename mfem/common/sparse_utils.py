try:
   import mfem.par as mfem
   MFEM_PAR = True
except:
   import mfem.ser as mfem    
   MFEM_PAR = False

from scipy.sparse import csr_matrix, coo_matrix, lil_matrix, csc_matrix

def eliminate_rows(m, rows):
    '''
    if m is not csr, it is converted to csr
    return csr_matrix
    '''
    if not isinstance(m, csr_matrix):
        try:
            m = m.tocsr()
        except:
            raise ValueError('Can not convert Matrix to CSR format.')            
    for row in rows:
        m.data[m.indptr[row]:m.indptr[row+1]] = 0.0  #  m[x,:] = 0.0
    m.eliminate_zeros()
    return m

def eliminate_cols(m, cols):
    '''
    if m is not csr, it is converted to csc
    return csc_matrix
    '''
    if not isinstance(m, csc_matrix):
        try:
            m = m.tocsc()
        except:
            raise ValueError('Can not convert Matrix to CSR format.')            
    for col in cols:
       m.data[m.indptr[col]:m.indptr[col+1]] = 0.0  #  m[:,x] = 0.0
    m.eliminate_zeros()
    return m

def sparsemat_to_scipycsr(mat, dtype):
     w, h = mat.Width(), mat.Height()
     I = mat.GetIArray()
     J = mat.GetJArray()
     data = mat.GetDataArray()
     mat.LoseData()
     m1 =csr_matrix((data, J, I), shape = (h, w),
                    dtype = dtype)
     return m1
        

