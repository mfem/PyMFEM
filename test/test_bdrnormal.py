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

    idx = mfem.intArray(list(np.unique(mesh.GetBdrAttributeArray())))

    gf.ProjectBdrCoefficient(func, idx)

    print("File: " + filename)
    print("-- maximum geometry error", np.max(gf.GetDataArray()))


def run_test():
    check_geometry_error('sphere-o2.mesh')
    check_geometry_error('sphere-o3.mesh')
    check_geometry_error('sphere-o4.mesh')
    check_geometry_error('sphere-o5.mesh')


if __name__ == '__main__':
    run_test()
