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


   example_4 tests
      how to make complex matrix
      transpose/conj of complex matrix

      mpirun -np 6 python2.7 example_4.py
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
shape = (1, 6)
# make sample matrix
row  = np.array([0, 0, ])
col  = np.array([0, 3, ])
data = np.array([1, 3.*(myid+1), ])
m = coo_matrix((data, (row, col)), shape=shape)
m = m.tocsr()


M1 = CHypreMat(m, m)
M2 = M1.transpose().conj()

print_hypre(M1[1], 'matrix M1 (imag)')
print_hypre(M2[1], 'matrix M2 (imag)')





