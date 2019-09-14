from __future__ import print_function
import os
from os.path import expanduser, join
import sys
import numpy as np

from mfem import path as mfem_path

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
    meshfile =expanduser(join(mfem_path, 'data', 'beam-tri.mesh'))
    mesh = mfem.Mesh(meshfile, 1,1)
    fec = mfem.H1_FECollection(1, 1)
    fespace = mfem.FiniteElementSpace(mesh, fec)

    fespace.Save()
    
if __name__=='__main__':
    run_test()
 
