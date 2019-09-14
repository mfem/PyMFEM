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
    a = mfem.intArray([1,2,3])

    for i in range(a.Size()):
        print(a[i])
    a.Print()
    a.Print("intArray.dat")
    a.Save("intArray_save.dat")
    
    b = mfem.doubleArray([1.1,2.2,3.])
    for i in range(b.Size()):
        print(b[i])
    b.Print()
    b.Print("doubleArray.dat")
    b.Save("doubleArray_save.dat")    
    
if __name__=='__main__':
    run_test()
 
