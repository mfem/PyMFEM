import sys
if 'mfem.ser' in sys.modules:
   import mfem.ser as mfem
else:
   import mfem.par as mfem
   
import numpy as np

def ToHypreParCSR(mat):
    '''
    convert scipy sparse matrix to hypre

    vertically stack csr matrix to generte HYPRE Par CSR
    '''
    from mpi4py import MPI

    if mfem.sizeof_HYPRE_Int() == 4:
        dtype = 'int32'
    else:
        dtype = 'int64'        
    

    mat = mat.astype('float')
    comm     = MPI.COMM_WORLD     
    num_proc = MPI.COMM_WORLD.size
    myid     = MPI.COMM_WORLD.rank

    ml, nl = mat.shape

    # collect row array
    m_array = comm.allgather(ml)
    n_array = comm.allgather(nl)    

    row_starts = [0] + list(np.cumsum(m_array))
    col_starts = [0] + list(np.cumsum(n_array))    
    row_starts = np.array([row_starts[myid], row_starts[myid+1], row_starts[-1]], dtype=dtype)
    col_starts = np.array([col_starts[0], col_starts[1], col_starts[1]], dtype=dtype)    
    m = row_starts[-1]
    n = col_starts[-1]
    nrows = ml

    i = mat.indptr.astype(dtype)
    j = mat.indices.astype(dtype)
    data = mat.data

    #
    # here I am passing row_starts as col_starts too
    # it seems row_starts and col_starts are both to determin
    # which part is treated diagnal element.
    #
    return  mfem.HypreParMatrix(MPI.COMM_WORLD,
                                nrows,
                                m, n, [i, j,
                                data, row_starts, row_starts])

def ToScipyCoo(mat):
    '''
    convert HypreParCSR to Scipy Coo Matrix
    '''
    num_rows, ilower, iupper, jlower, jupper, irn, jcn, data = mat.GetCooDataArray()
    from mpi4py import MPI
    myid     = MPI.COMM_WORLD.rank
    
    m = iupper - ilower + 1
    n = jupper - jlower + 1
    n = mat.N()    

    #if myid == 1:
    #   print m, n
    #   print data
    #   print irn
    #   print jcn
    from scipy.sparse import coo_matrix
    return coo_matrix((data, (irn-ilower, jcn)), shape = (m, n))

def ParAdd(A, B):
    '''
    add HypreParCSR
    '''
    return ToHypreParCSR((ToScipyCoo(A)+ ToScipyCoo(B)).tocsr())
    
def ParMultComplex(A, B):
    '''
    compute complex mult of hypre real matrices

    A = (R_A, I_A)
    B = (R_B, I_B)

    (R_A*R_B - I_A*I_B, R_A*I_B + I_A*R_B)
    '''
    R_A, I_A = A
    R_B, I_B = B

    if I_A is None and I_B is None:
       return (mfem.ParMult(R_A, R_B), None)
    elif I_A is None:
       r = mfem.ParMult(R_A, R_B)
       i = mfem.ParMult(R_A, I_B)
       return (r, i)
    elif I_B is None:
       r = mfem.ParMult(R_A, R_B)
       i = mfem.ParMult(I_A, R_B)
       return (r, i)       
    else:
       r = ToHypreParCSR((ToScipyCoo(mfem.ParMult(R_A, R_B)) - ToScipyCoo(mfem.ParMult((I_A, I_B)))).tocsr())
       i = ToHypreParCSR((ToScipyCoo(mfem.ParMult(R_A, I_B)) + ToScipyCoo(mfem.ParMult((I_A, R_B)))).tocsr())
       return (r, i)

def TransposeComplex(A):
    '''
    A is tuple (A_real, A_imag), whose real/imag are
    HypreParCSR
    '''
    return (A[0].Transpose(), A[1].Transpose())
 
def RapComplex(A, B):
    '''
    Bt * A * B

    for complex A and B
    '''
    return ParMultComplex(TransposeComplex(B), ParMultComplex(A, B))
