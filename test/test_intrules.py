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
    a = mfem.IntegrationPointArray()
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.IntegrationPointArray(mt)
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.IntegrationPointArray(3, mt)
    assert a.Size() == 3, "array size is wrong"
    
    ir = mfem.IntegrationRule(5)

    pts = [ir[i] for i in range(5)]
    a = mfem.IntegrationPointArray(pts)

    for i in range(5):
        print(a[i])

    try:
        a = mfem.IntegrationPointArray([1,2,3])
    except:
        print("this one is supposed to fail")
        
if __name__=='__main__':
    run_test()
