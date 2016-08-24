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

    print mat.shape
    print m_array
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

    return  mfem.HypreParMatrix(MPI.COMM_WORLD,
                                nrows,
                                m, n, [i, j,
                                       data, row_starts, col_starts])

