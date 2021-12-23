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
    a = mfem.GeometryTypeArray()
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.GeometryTypeArray(mt)
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.GeometryTypeArray(3, mt)
    assert a.Size() == 3, "array size is wrong"
    
    mesh = mfem.Mesh(5, 5, "TRIANGLE")


    gg1 = [mesh.GetElementGeometry(i) for i in range(mesh.GetNE())]
    a = mfem.GeometryTypeArray(gg1)
    gg2 = [a[i] for i in range(a.Size())]

    assert gg1 == gg2, "does not match"
    
if __name__=='__main__':
    run_test()
