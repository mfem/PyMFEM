from __future__ import print_function
import os
import sys
import numpy as np

from scipy.sparse import csr_matrix, coo_matrix

def run_test():
    import mfem.par as par
    from mfem.common.parcsr_extra import ToHypreParCSR, ToScipyCoo
    from mpi4py import MPI
    
    comm     = MPI.COMM_WORLD     
    num_proc = MPI.COMM_WORLD.size
    myid     = MPI.COMM_WORLD.rank

    def print_hypre(M, txt):
        for i in range(num_proc):
            MPI.COMM_WORLD.Barrier()                              
            if myid == i:
                if myid == 0:
                    print(txt)
                    print('MyID: ', myid)
                else:
                    print('MyID: ', myid)
                print(ToScipyCoo(M))

    # make sample matrix
    row  = np.array([0, 0, 1, 1])
    col  = np.array([0, 3, 1, 2])
    data = np.array([4, 5, 7, 9])
    m = coo_matrix((data, (row, col)), shape=(2,4))
    m = m.tocsr()
    m = m*(myid+1)

    M = ToHypreParCSR(m, assert_non_square_no_col_starts=False)
    print_hypre(M, 'matrix M')
    
if __name__=='__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':   
        import mfem.par as mfem
        run_test()        
    else:
        import mfem.ser as mfem
    

