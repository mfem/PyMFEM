import sys
import numpy as np
from scipy.sparse import csr_matrix

if len(sys.argv) < 2:
    import mfem.ser as mfem
elif sys.argv[1] == 'par':
    import mfem.par as mfem    


smat = csr_matrix([[1,2,3,], [0, 0, 1], [2, 0, 0]])
mmat = mfem.SparseMatrix(smat)

offset1 = mfem.intArray([0, 3, 6])
offset2 = mfem.intArray([0, 3, 6])

m = mfem.BlockMatrix(offset1, offset2)
m.SetBlock(0, 0, mmat)

print(m._offsets)
print(m._linked_mat)

m = mfem.BlockMatrix(offset1)
m.SetBlock(0, 1, mmat)
m.SetBlock(1, 0, mmat)
print(m._offsets)
print(m._linked_mat)
