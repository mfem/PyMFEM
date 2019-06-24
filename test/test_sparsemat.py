from __future__ import print_function
import os
import sys
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == '-p':
   print("loading parallel")
   import mfem.par as mfem
else:
   import mfem.ser as mfem
        

def run_test():
    print("Test sparsemat module")
    indptr = np.array([0, 2, 3, 6], dtype=np.int32)
    indices = np.array([0, 2, 2, 0, 1, 2], dtype=np.int32)
    data = np.array([1, 2, 3, 4, 5, 6], dtype=float)
    
    mat = mfem.SparseMatrix([indptr, indices, data, 3,3])
    mat.Print()
    mat.PrintInfo()
    
if __name__=='__main__':
    run_test()
