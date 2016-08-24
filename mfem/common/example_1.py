import mfem.par as par
from mfem.common.scipycsr_to_hypreparcsr import ToHypreParCSR


from scipy.sparse import csr_matrix

from mpi4py import MPI


comm     = MPI.COMM_WORLD     
num_proc = MPI.COMM_WORLD.size
myid     = MPI.COMM_WORLD.rank

m = csr_matrix((4,2))
m[myid, myid] = 1.0

M = ToHypreParCSR(m)

M.Print('matrix')

