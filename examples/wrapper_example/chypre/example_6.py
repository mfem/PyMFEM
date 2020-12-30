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

row  = np.array([0, 0, ])
col  = np.array([0, myid+1, ])
data = np.array([1, myid+1, ])
m = coo_matrix((data, (row, col)), shape=shape)
m = m.tocsr()

M1 = CHypreMat(m, m*2)

print_hypre(M1.real,      '#### matrix M1 real ')
print_hypre(M1.imag,      '#### matrix M1 imag ')
print_hypre((-M1).real,   '#### matrix -M1 real')
print_hypre((-M1).imag,   '#### matrix -M1 imag')
print_hypre((M1-M1).real, '#### M1-M2 (real)')
print_hypre((M1-M1).imag, '#### M1-M2 (imag)')

print_hypre(M1.transpose().real, '#### matrix M1^t real ')






