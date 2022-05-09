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
    a = mfem.intArray()
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.intArray(mt)
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.intArray(3, mt)
    assert a.Size() == 3, "array size is wrong"

    
    a = mfem.intArray([1,2,3])
    assert a.Size() == 3, "array size is wrong"
    for i in range(a.Size()):
        assert a[i] == i + 1, "wrong data is set"

    a.Print()
    a.Print("intArray.dat")
    a.Save("intArray_save.dat")

    a2 = mfem.intArray([a.GetData(), a.Size()])
    a3 = mfem.intArray((a.GetData(), a.Size()))    

    a.ToList()
    a[2] = 4
    assert a2.ToList() == a3.ToList(), "pointer constructer fails"
    
    b = mfem.doubleArray([1.1,2.2,3.])
    for i in range(b.Size()):
        print(b[i])
    b.Print()
    b.Print("doubleArray.dat")
    b.Save("doubleArray_save.dat")

    b2 = mfem.doubleArray([b.GetData(), b.Size()])
    b3 = mfem.doubleArray((b.GetData(), b.Size()))
    b[1] = 4    
    print(b2.ToList())
    print(b3.ToList())    

    a = mfem.boolArray([True]*3)
    print(a[1])
    print(a.ToList())    
    assert a.Size() == 3, "array size is wrong"
    
    
if __name__=='__main__':
    run_test()
 
