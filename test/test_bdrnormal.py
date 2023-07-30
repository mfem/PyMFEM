'''
 test BondaryNormalCoefficient

'''
from numba import cfunc, farray, carray
import os
from os.path import expanduser, join
import sys
import io
from mfem import path as mfem_path
import numpy as np
import numba
import time
from math import sqrt


if len(sys.argv) > 1 and sys.argv[1] == '-p':
    import mfem.par as mfem
    use_parallel = True
    from mfem.common.mpi_debug import nicePrint
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0
    nicePrint = print


def check_geometry_error(filename):
    meshfile = expanduser(join("..", 'data', filename))
    mesh = mfem.Mesh(meshfile)

    print("File: " + filename)
    print("NE", mesh.GetNE())
    print("NBE", mesh.GetNBE())
    norm = mfem.VectorBdrNormalCoefficient(mesh.SpaceDimension())

    @mfem.jit.scalar(dependency=(norm,))
    def func(ptx, norm):
        #print("point", ptx, norm)
        #print("error", np.sum(ptx - norm)**2, sqrt(np.sum(ptx**2)))
        #print(sqrt(np.sum(ptx**2)), np.sum(norm**2))
        return sqrt(np.sum(ptx - norm)**2)

    order = 1
    dim = mesh.Dimension()
    fec = mfem.H1_FECollection(order, dim)
    fes = mfem.FiniteElementSpace(mesh, fec, 1)

    gf = mfem.GridFunction(fes)
    gf.Assign(0.0)


    if mesh.Dimension() == 3:
        idx = mfem.intArray(list(np.unique(mesh.GetBdrAttributeArray())))
        gf.ProjectBdrCoefficient(func, idx)
    else:
        gf.ProjectCoefficient(func)

    print("-- maximum geometry error", np.max(gf.GetDataArray()))

    gfname = filename[:-5]+'-data.gf'
    gf.Save(gfname)
    print("-- error is saved as " + gfname)

def run_test():
    check_geometry_error('sphere-o2.mesh')
    check_geometry_error('sphere-o3.mesh')
    check_geometry_error('sphere-o4.mesh')
    check_geometry_error('sphere-o5.mesh')
    check_geometry_error('sphere-surface-o3.mesh')


if __name__ == '__main__':
    run_test()
