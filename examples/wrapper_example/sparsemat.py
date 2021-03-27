import numpy as np
from scipy.sparse import csr_matrix

import mfem.ser as mfem


smat = csr_matrix([[1,2,3,], [0, 0, 1], [2, 0, 0]])
print(mfem.SparseMatrix(smat))
print(smat)
mmat = mfem.SparseMatrix(smat)
mmat.Print()
