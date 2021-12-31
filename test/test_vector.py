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
    a = mfem.Vector(np.arange(20.0))
    a.Print()
    
    a *= 30.
    a.Print()

    a[-15:].Print()
    a[-15:].Print_HYPRE()    

    a.Print("vector.dat")
    a.Print_HYPRE("vector_hypre.dat")


    x = mfem.VectorPtrArray([a]*3)
    x[2].Print()
    
if __name__=='__main__':
    run_test()
 
