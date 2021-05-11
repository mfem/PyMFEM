from __future__ import print_function
import os
from os.path import expanduser, join
import sys
import io
from mfem import path as mfem_path
import numba

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

from numba import cfunc, carray    

@cfunc(mfem.scalar_sig)
def s_func(ptx, sdim):
    return 2

def run_test():
    #meshfile = expanduser(join(mfem_path, 'data', 'semi_circle.mesh'))
    mesh = mfem.Mesh(6, 6, "TRIANGLE")
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()    
    fec = mfem.H1_FECollection(1, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec, 1)
    print('Number of finite element unknowns: '+
          str(fespace.GetTrueVSize()))

    c = mfem.NumbaFunction(s_func, sdim).GenerateCoefficient()    
    gf = mfem.GridFunction(fespace)
    gf.ProjectCoefficient(c)

    gf.Print()
if __name__=='__main__':
    run_test()
    
