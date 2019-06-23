from __future__ import print_function
import os
import sys
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == '-p':   
   import mfem.par as mfem
else:
   import mfem.ser as mfem
        

def run_test():
    print("Test sparsemat module")
    indptr = np.array([0, 2, 3, 6], dtype=int)
    indices = np.array([0, 2, 2, 0, 1, 2])
    data = np.array([1, 2, 3, 4, 5, 6], dtype=float)
    
    mat = mfem.SparseMatrix([indptr, indices, data, 2,2])
    mat.Print()

    print(mat.GetIArray())
    print(mat.GetJArray())
    print(mat.GetDataArray())    
        
    
if __name__=='__main__':
    run_test()
