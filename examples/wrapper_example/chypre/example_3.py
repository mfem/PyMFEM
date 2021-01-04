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


   example_3 tests
      multiplication between square matrix and non-square matrix.
   
      it seems,,
      1) m>=n is requred to generate ParCSR ??
      2) partitioning of diagnal elemnts should be the same,
         which makes sense since diagnals are computed w/o
         communication.



'''
import mfem.par as par
from mfem.common.parcsr_extra import *

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
shape = (2, 5) if myid < 2 else (1, 5)
# make sample matrix
row  = np.array([0, 0, 0, 0])
col  = np.array([0, 3, 1, 2])
data = np.array([4, 5, 7, 9])
m = coo_matrix((data, (row, col)), shape=shape)
m = m.tocsr()
m = m*(myid+1)

M = ToHypreParCSR(m)
print_hypre(M, 'matrix M')

# make sample vector, it produce vertical array
M2 = Array2Hypre(np.array([1,2,3,4,5]))
'''
if myid == 0:
   row  = np.array([0, 1, ])
   col  = np.array([0, 0, ])
   data = np.array([4, 10,  ])
   m2 = coo_matrix((data, (row, col)), shape=(2, 3))
   m2 = m2.tocsr()
elif myid == 1:
   row  = np.array([ ])
   col  = np.array([ ])
   data = np.array([ ])
   m2 = coo_matrix((data, (row, col)), shape=(2, 3))
   m2 = m2.tocsr()
#   m2 = None
else:
   row  = np.array([])
   col  = np.array([])
   data = np.array([])
   m2 = coo_matrix((data, (row, col)), shape=(1, 3))
   m2 = m2.tocsr()
M2 = ToHypreParCSR(m2)
'''
#print ToScipyCoo(M2)

# adding matrix
M3 = mfem.ParMult(M, M2)
print_hypre(M3, 'M*v')
v = Hypre2Array(M3)
if myid == 0: print(v)
M3 = mfem.ParMult(M2.Transpose(), M)
print_hypre(M3, 'v^t * M')
v = Hypre2Array(M3)
if myid == 0: print(v)
