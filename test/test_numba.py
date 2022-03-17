from __future__ import print_function
from numba import cfunc, farray, carray
import os
from os.path import expanduser, join
import sys
import io
from mfem import path as mfem_path
import numpy as np
import numba
import time


if len(sys.argv) > 1 and sys.argv[1] == '-p':
    import mfem.par as mfem
    use_parallel = True
    from mfem.common.mpi_debug import nicePrint as print
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0


class s_coeff(mfem.PyCoefficient):
    def __init__(self):
        mfem.PyCoefficient.__init__(self)

    def EvalValue(self, p):
        return p[0]


class v_coeff(mfem.VectorPyCoefficient):
    def __init__(self, dim):
        mfem.VectorPyCoefficient.__init__(self, dim)

    def EvalValue(self, p):
        return (p[0], p[1], p[2])


class m_coeff(mfem.MatrixPyCoefficient):
    def __init__(self, dim):
        mfem.MatrixPyCoefficient.__init__(self, dim)

    def EvalValue(self, p):
        return np.array([[p[0], p[1], p[2]],
                         [0.0, p[1], p[2]],
                         [0.0, 0.0, p[2]]])

@cfunc("float64(float64, float64, float64)")
def s_func0(x, y, z):
    return x

def s_func1(x, y, z):
    return x
s_func0 = cfunc("float64(float64, float64, float64)")(s_func1)

@cfunc(mfem.scalar_sig, cache=False)
def s_func(ptx, sdim):
    return s_func0(ptx[0], ptx[1],  ptx[2])


@cfunc(mfem.vector_sig)
def v_func(ptx, out, sdim, vdim):
    out_array = carray(out, (vdim, ))
    for i in range(sdim):
        out_array[i] = ptx[i]


@cfunc(mfem.matrix_sig)
def m_func(ptx, out, sdim, vdim):
    # we use farray to assign the data like above.
    # note out is zero-ed in wrapper. so we don't need to
    # set zero here.
    out_array = farray(out, (vdim, vdim))
    out_array[0, 0] = ptx[0]
    out_array[0, 1] = ptx[1]
    out_array[0, 2] = ptx[2]
    #out_array[1, 0] = 0.0
    out_array[1, 1] = ptx[1]
    out_array[1, 2] = ptx[2]
    #out_array[2, 0] = 0.0
    #out_array[2, 1] = 0.0
    out_array[2, 2] = ptx[2]
    '''
    accessing the array linearly does not speed up.
    out_array = carray(out, (vdim * vdim, ))
    out_array[0] = ptx[0]
    out_array[1] = 0.0
    out_array[2] = 0.0
    out_array[3] = ptx[1]
    out_array[4] = ptx[1]
    out_array[5] = 0.0
    out_array[6] = ptx[2]
    out_array[7] = ptx[2]
    out_array[8] = ptx[2]
    '''

# if dim is know, this provides a simpler way to use matrix coefficient
@mfem.jit.matrix()
def m_func2(ptx, out):
    out_array = farray(out, (3, 3))
    out_array[0, 0] = ptx[0]
    out_array[0, 1] = ptx[1]
    out_array[0, 2] = ptx[2]
    out_array[1, 1] = ptx[1]
    out_array[1, 2] = ptx[2]
    out_array[2, 2] = ptx[2]
    
def check(a, b, msg):
    assert len(a) == len(b), msg
    assert np.sum(np.abs(a - b)) == 0, msg


def run_test():
    #meshfile = expanduser(join(mfem_path, 'data', 'semi_circle.mesh'))
    mesh = mfem.Mesh(3, 3, 3, "TETRAHEDRON")
    mesh.ReorientTetMesh()

    order = 1

    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()
    fec1 = mfem.H1_FECollection(order, dim)
    fespace1 = mfem.FiniteElementSpace(mesh, fec1, 1)

    fec2 = mfem.ND_FECollection(order, dim)
    fespace2 = mfem.FiniteElementSpace(mesh, fec2, 1)

    print("Element order :", order)

    print('Number of H1 finite element unknowns: ' +
          str(fespace1.GetTrueVSize()))
    print('Number of ND finite element unknowns: ' +
          str(fespace2.GetTrueVSize()))

    print("Checking scalar")

    gf = mfem.GridFunction(fespace1)
    c1 = mfem.NumbaFunction(s_func, sdim).GenerateCoefficient()

    @mfem.jit.scalar(sdim)
    def c11(ptx, _sdim):
        return s_func0(ptx[0], ptx[1],  ptx[2])
    @mfem.jit.scalar()
    def c12(ptx):
        return s_func0(ptx[0], ptx[1],  ptx[2])

    c2 = s_coeff()

    gf.Assign(0.0)
    start = time.time()
    gf.ProjectCoefficient(c12)
    end = time.time()
    data1 = gf.GetDataArray().copy()
    print("Numba time (scalar)", end - start)

    gf.Assign(0.0)
    start = time.time()
    gf.ProjectCoefficient(c2)
    end = time.time()
    data2 = gf.GetDataArray().copy()
    print("Python time (scalar)", end - start)

    check(data1, data2, "scalar coefficient does not agree with original")

    print("Checking vector")
    gf = mfem.GridFunction(fespace2)
    c3 = mfem.VectorNumbaFunction(v_func, sdim, dim).GenerateCoefficient()
    c4 = v_coeff(dim)

    gf.Assign(0.0)
    start = time.time()
    gf.ProjectCoefficient(c3)
    end = time.time()
    data1 = gf.GetDataArray().copy()
    print("Numba time (vector)", end - start)

    gf.Assign(0.0)
    start = time.time()
    gf.ProjectCoefficient(c4)
    end = time.time()
    data2 = gf.GetDataArray().copy()
    print("Python time (vector)", end - start)

    check(data1, data2, "vector coefficient does not agree with original")

    print("Checking matrix")
    a1 = mfem.BilinearForm(fespace2)
    a2 = mfem.BilinearForm(fespace2)
    a3 = mfem.BilinearForm(fespace2)    
    c4 = mfem.MatrixNumbaFunction(m_func, sdim, dim).GenerateCoefficient()
    c5 = m_coeff(dim)

    a1.AddDomainIntegrator(mfem.VectorFEMassIntegrator(c4))
    a2.AddDomainIntegrator(mfem.VectorFEMassIntegrator(c5))
    a3.AddDomainIntegrator(mfem.VectorFEMassIntegrator(m_func2))
    
    start = time.time()
    a1.Assemble()
    end = time.time()
    a1.Finalize()
    M1 = a1.SpMat()

    print("Numba time (matrix)", end - start)

    start = time.time()
    a2.Assemble()
    end = time.time()
    a2.Finalize()
    M2 = a2.SpMat()
    print("Python time (matrix)", end - start)

    start = time.time()
    a3.Assemble()
    end = time.time()
    a3.Finalize()
    M3 = a3.SpMat()
    print("Numba (simpler interface) (matrix)", end - start)
    

    #from mfem.commmon.sparse_utils import sparsemat_to_scipycsr
    #csr1 = sparsemat_to_scipycsr(M1, float)
    #csr2 = sparsemat_to_scipycsr(M2, float)

    check(M1.GetDataArray(),
          M2.GetDataArray(),
          "matrix coefficient does not agree with original")
    check(M1.GetIArray(),
          M2.GetIArray(),
          "matrix coefficient does not agree with original")
    check(M1.GetJArray(),
          M2.GetJArray(),
          "matrix coefficient does not agree with original")
    check(M1.GetDataArray(),
          M3.GetDataArray(),
          "matrix coefficient does not agree with original")
    check(M1.GetIArray(),
          M3.GetIArray(),
          "matrix coefficient does not agree with original")
    check(M1.GetJArray(),
          M3.GetJArray(),
          "matrix coefficient does not agree with original")

    print("PASSED")


if __name__ == '__main__':
    run_test()
