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
    a = mfem.uintArray(5)
    b = mfem.int8Array(5)
    c = mfem.int64Array(5)
    d = mfem.boolArray(5)

    def check(a, value):
        a.GetDataArray()[:] = value
        print(a[0], a[-1], type(a[0]), )
        assert a[0] == value, "Array data element does not agree"


    check(a, 5)
    check(b, 100)
    check(c, 1000)
    check(d, True)    
    
    
if __name__=='__main__':
    run_test()
