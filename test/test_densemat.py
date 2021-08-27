from __future__ import print_function
import os
import sys
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == '-p':   
    import mfem.par as mfem
    use_parallel = True
    from mfem.common.mpi_debug import nicePrint as print
    from mpi4py import MPI
    myid  = MPI.COMM_WORLD.rank
    
else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0
    
def run_test(): 
    m = mfem.DenseMatrix(3,3)
    m.Assign(np.arange(9.).reshape(3,3))
    m.Print()
    print(np.arange(9.).reshape(3,3))

    m = mfem.DenseMatrix(3,5)
    m.Assign(np.arange(15.).reshape(3,5))
    m.Print()
    print(np.arange(15.).reshape(3,5))
    
    
    x = np.zeros(5)+1
    y = np.zeros(5)+2
    z = np.zeros(5)+3
    mm = mfem.DenseMatrix(3, 5)
    mm.Assign(np.vstack([x, y, z]))
    mm.Print()
    mfem.Vector(mm.GetData(), 15).Print()

    mm.Print("densmat.dat")

    ### added this check based on the issue 99.
    import tracemalloc
    import gc
    tracemalloc.start()
    n = 4 * 1024 // 8
    d = 1024
    K = mfem.DenseMatrix(d,d)
    for i in range(10): 
        K.Assign(np.diag(d * np.ones(d)))
        data = K.GetDataArray()
        size, peak = tracemalloc.get_traced_memory()        
        #print(f"{size=}, {peak=}")

    m = mfem.DenseTensor(2, 3, 5)
    m.Assign(0)
    m[1, 2, 3] = 5
    print(m[1, 2, 3])
    
    assert np.all(m[3].GetDataArray() == m.GetDataArray()[3]), "First column slicing does not agree"

    t = np.arange(30.).reshape(2, 3, 5)
    size, peak = tracemalloc.get_traced_memory()        
    print(f"{size=}, {peak=}")
    
    m.Assign(t)
    data = []
    for i in range(m.SizeI()):
       for j in range(m.SizeJ()):
        for k in range(m.SizeK()):
            assert m[i, j, k] == t[i, j, k], "Tensor assigne fails"

    # these two should not agree (internal memory layout differs)        
    print(t.flatten())
    print(mfem.Vector(m.Data(), 30).GetDataArray())
    
if __name__=='__main__':
    run_test()
 
