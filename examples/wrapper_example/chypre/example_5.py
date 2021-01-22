'''
 This example shows
   constructing HYPREParCSR in a distributed manner
     1) build scipy.sparse on each node
           each node has informaiton of small number of rows
     2) assemble it to HYPREParCSR
           scipy matries are concatenated vertical

   convert back ParCSR to scipy.sparse
     1) on each node, it produce scipy.sparse which contain
        only its own rows.


   example_5 tests
      how to make a non-square matrix
      transpose/conj of complex matrix

   usage: mpirun -np 3 python example_5.py

'''
import mfem.par as par
from mfem.common.chypre import *

from scipy.sparse import csr_matrix, coo_matrix

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
          m  = ToScipyCoo(M)              
          print('shape = ', m.shape)
          print(m.toarray())
       MPI.COMM_WORLD.Barrier()                              
shape = (2, 10)
# make sample matrix
col_starts = [0, 3, 7, 10]

row  = np.array([0, 0, ])
col  = np.array([0, 3*(myid+1), ])
data = np.array([1, 3.*(myid+1), ])
m = coo_matrix((data, (row, col)), shape=shape)
m = m.tocsr()

M1 = CHypreMat(m, None, col_starts = 
       [col_starts[myid], col_starts[myid+1], col_starts[-1]])
M2 = M1.transpose().conj()

print_hypre(M1[0], '#### matrix M1 ')
print_hypre(M2[0], '#### matrix M2 (transposed)')

print_hypre((M1*M2)[0], '#### matrix M1*M2 (6x10)*(10x6) => 6x6')
print_hypre((M2*M1)[0], '#### matrix M2*M1 (10x6)*(6x10) => 10x10')



