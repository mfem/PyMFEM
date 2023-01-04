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


@mfem.jit.matrix(sdim=3)
def m_func2(p,):
    #out_array = farray(out, (3, 3))
    #out_array[0, 0] = ptx[0]
    #out_array[0, 1] = ptx[1]
    #out_array[0, 2] = ptx[2]
    #out_array[1, 1] = ptx[1]
    #out_array[1, 2] = ptx[2]
    #out_array[2, 2] = ptx[2]
    return np.array([[p[0], p[1], p[2]],
                     [0.0, p[1], p[2]],
                     [0.0, 0.0, p[2]]])


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
    def c11(ptx):
        return s_func0(ptx[0], ptx[1],  ptx[2])

    @mfem.jit.scalar(td=True, debug=True)
    def c12(ptx, t):
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

    @mfem.jit.vector(sdim=3, dependency=(c3, c11), td=True,  complex=True)
    def v_func4(ptx, t, c3, c11):
        return np.array([c3[0], c3[1], c3[2]], dtype=np.complex128)
    # @mfem.jit.vector(sdim=3, complex=True)
    # def v_func4(ptx,  ):
    #    return np.array([ptx[0],ptx[1],ptx[2]], dtype=np.complex128)

    gf.Assign(0.0)
    start = time.time()
    for i in range(10):
        gf.ProjectCoefficient(c3)
    end = time.time()
    data1 = gf.GetDataArray().copy()
    print("Numba time (vector)", end - start)

    gf.Assign(0.0)
    start = time.time()
    for i in range(10):
        gf.ProjectCoefficient(c4)
    end = time.time()
    data2 = gf.GetDataArray().copy()
    print("Python time (vector)", end - start)

    gf.Assign(0.0)
    start = time.time()
    for i in range(10):
        gf.ProjectCoefficient(v_func4.real)
    end = time.time()
    data3 = gf.GetDataArray().copy()
    print("Numba2 time (vector)", end - start)

    check(data1, data2, "vector coefficient does not agree with original")
    check(data1, data3, "vector coefficient does not agree with original")

    print("speed comparision with C++")

    @mfem.jit.vector(sdim=3, interface="c++", debug=True)
    def v_func4_old(ptx, out):
        out[0] = 1
        out[1] = 2.
        out[2] = 3

    @mfem.jit.vector(sdim=3)
    def v_func4_new(ptx):
        return np.array([1, 2, 3.])

    gf.Assign(0.0)
    start = time.time()
    for i in range(10):
        gf.ProjectCoefficient(v_func4_old)
    end = time.time()
    data3 = gf.GetDataArray().copy()
    print("Numba time (vector) - old interface", end - start)

    gf.Assign(0.0)
    start = time.time()
    for i in range(10):
        gf.ProjectCoefficient(v_func4_new)
    end = time.time()
    data3 = gf.GetDataArray().copy()
    print("Numba time (vector) - new interface", end - start)

    val = mfem.Vector([1, 2, 3])
    cc = mfem.VectorConstantCoefficient(val)
    gf.Assign(0.0)
    start = time.time()
    for i in range(10):
        gf.ProjectCoefficient(cc)
    end = time.time()
    data3 = gf.GetDataArray().copy()
    print("C++ constant coefficient", end - start)

    print("Checking matrix")
    a1 = mfem.BilinearForm(fespace2)
    a2 = mfem.BilinearForm(fespace2)
    a3 = mfem.BilinearForm(fespace2)
    a4 = mfem.BilinearForm(fespace2)
    a5 = mfem.BilinearForm(fespace2)

    c4 = mfem.MatrixNumbaFunction(m_func, sdim, dim).GenerateCoefficient()
    c5 = m_coeff(dim)

    a1.AddDomainIntegrator(mfem.VectorFEMassIntegrator(c4))
    a2.AddDomainIntegrator(mfem.VectorFEMassIntegrator(c5))
    a3.AddDomainIntegrator(mfem.VectorFEMassIntegrator(m_func2))

    @mfem.jit.matrix(sdim=3, dependency=(c4, c4), complex=True, td=True)
    def m_func3(ptx, t, c4, c5):
        ret = np.array([[c5[0, 0], c5[0, 1], c5[0, 2]],
                        [0.0,      c5[1, 1], c5[1, 2]],
                        [0.0,      0.0,      c5[2, 2]]])
        return (t/3.0)*ret*1j

    @mfem.jit.matrix(sdim=3, dependency=(m_func3, c4))
    def m_func4_complex(ptx, m_func3, c5):
        return m_func3.imag

    @mfem.jit.matrix(sdim=3, dependency=((m_func3.real, m_func3.imag), c4), td=True)
    def m_func4_split(ptx, t, m_func3, c5):
        return m_func3.imag*(t/3.0)

    '''
    @mfem.jit.matrix(sdim=3, complex=False, debug=True)
    def m_func5(p):
        x = p[0]
        y = p[1]
        return np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    '''
    # _out_ =
    # return np.array(_out_)

    a4.AddDomainIntegrator(mfem.VectorFEMassIntegrator(m_func4_complex))
    a5.AddDomainIntegrator(mfem.VectorFEMassIntegrator(m_func4_split))

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

    start = time.time()
    m_func4_complex.SetTime(3.0)
    a4.Assemble()
    end = time.time()
    a4.Finalize()
    M4 = a4.SpMat()
    print("Numba (complex dependency as complex) (matrix)", end - start)

    start = time.time()
    m_func4_split.SetTime(3.0)
    a5.Assemble()
    end = time.time()
    a5.Finalize()
    M5 = a5.SpMat()
    print("Numba (complex dependency as decomposed) (matrix)", end - start)

    #from mfem.commmon.sparse_utils import sparsemat_to_scipycsr
    #csr1 = sparsemat_to_scipycsr(M1, float)
    #csr2 = sparsemat_to_scipycsr(M2, float)

    def compare_mat(M1o, M2o):
        check(M1o.GetDataArray(),
              M2o.GetDataArray(),
              "matrix coefficient does not agree with original")
        check(M1o.GetIArray(),
              M2o.GetIArray(),
              "matrix coefficient does not agree with original")
        check(M1o.GetJArray(),
              M2o.GetJArray(),
              "matrix coefficient does not agree with original")

    compare_mat(M1, M2)
    compare_mat(M1, M3)
    compare_mat(M1, M4)
    compare_mat(M1, M5)

    # print(m_func3.SpaceDimension())
    print("PASSED")


if __name__ == '__main__':
    run_test()
