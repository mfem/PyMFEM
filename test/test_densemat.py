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
    
    x = np.zeros(5)+1
    y = np.zeros(5)+2
    z = np.zeros(5)+3
    mm = mfem.DenseMatrix(3, 5)
    mm.Assign(np.vstack([x, y, z]))
    mm.Print()
    mfem.Vector(mm.GetData(), 15).Print()

    mm.Print("densmat.dat")
    
if __name__=='__main__':
    run_test()
 
