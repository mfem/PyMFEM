"""Native operators and Add/Subtract/Mult helpers; use -p for parallel bindings."""

import sys

import numpy as np
from numba import njit
from numba.core.errors import TypingError
import nsb_rt as nsb

if "-p" in sys.argv:
    import mfem.par as mfem
else:
    import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def vector_operations(a, b, output):
    a[0] = 3.0
    original = a[0]
    a += b
    a -= b
    a += b
    a -= 1.0
    a *= 2.0
    a /= 2.0
    a *= b
    a /= b
    mfem.Add(a, b, output)
    assert output[0] == 6.0
    mfem.Add(a, 2.0, b, output)
    assert output[0] == 8.0
    mfem.Add(2.0, a, b, output)
    assert output[0] == 12.0
    mfem.Add(2.0, a, 3.0, b, output)
    return original, a * b


@njit
def vector_subtraction(a, b, output):
    first = mfem.Subtract(a, b, output)
    difference = output.GetDataArray().copy()
    second = mfem.Subtract(2.0, a, b, output)
    return first, second, difference


@njit
def matrix_operations(a, b, output):
    inner = a * b
    a += b
    a -= b
    a *= 2.0
    mfem.Add(a, b, 3.0, output)
    assert output.GetDataArray()[0, 0] == 8.0
    mfem.Add(2.0, a, 3.0, b, output)
    assert output.GetDataArray()[0, 0] == 10.0
    mfem.Mult(a, b, output)
    return inner


@njit
def sparse_add(a, b):
    result = nsb.require_not_none(mfem.Add(a, b))
    result += a
    result *= 2.0
    return nsb.own_this(result)


@njit
def sparse_scaled_add(a, b):
    return nsb.own_this(nsb.require_not_none(mfem.Add(2.0, a, 3.0, b)))


@njit
def sparse_dense_add(a, output):
    mfem.Add(a, 2.0, output)


@njit
def sparse_dense_product(a, b):
    return nsb.own_this(nsb.require_not_none(mfem.Mult(a, b)))


@njit
def block_product(a, b):
    return nsb.own_this(nsb.require_not_none(mfem.Mult(a, b)))


@njit
def reflected(a):
    return 2.0 * a


@njit
def deferred(a, b):
    return mfem.Mult(a, b)


@njit
def deferred_none(a, b):
    return mfem.Mult(a, b, None)


@njit
def deferred_output(a, b, output):
    return mfem.Mult(a, b, output)


@njit
def hypre_operations(a, b):
    result = nsb.require_not_none(mfem.Add(2.0, a, 3.0, b))
    result += a
    result *= 2.0
    return nsb.own_this(result)


def run_test():
    assert not hasattr(mfem, "add_vector")
    assert not hasattr(mfem, "subtract_vector")
    for subtract in (vector_subtraction.py_func, vector_subtraction):
        x = mfem.Vector([3.0, 5.0, 9.0])
        y = mfem.Vector([1.0, 2.0, 4.0])
        difference = mfem.Vector(3)
        first, second, unscaled = subtract(x, y, difference)
        assert first is None and second is None
        np.testing.assert_allclose(unscaled, [2.0, 3.0, 5.0])
        np.testing.assert_allclose(difference.GetDataArray(), [4.0, 6.0, 10.0])
        np.testing.assert_allclose(x.GetDataArray(), [3.0, 5.0, 9.0])
        np.testing.assert_allclose(y.GetDataArray(), [1.0, 2.0, 4.0])
    a, b, output = mfem.Vector(3), mfem.Vector(3), mfem.Vector(3)
    a.Assign(2.0)
    b.Assign(2.0)
    assert vector_operations(a, b, output) == (3.0, 20.0)
    np.testing.assert_allclose(a.GetDataArray(), [4.0, 3.0, 3.0])
    np.testing.assert_allclose(output.GetDataArray(), [14.0, 12.0, 12.0])
    for unsupported, args in [(reflected, (a,))]:
        try:
            unsupported(*args)
        except TypingError:
            pass
        else:
            raise AssertionError("unsupported reflected operator compiled")

    da, db, dc = mfem.DenseMatrix(2), mfem.DenseMatrix(2), mfem.DenseMatrix(2)
    da.Assign(1.0)
    db.Assign(2.0)
    assert matrix_operations(da, db, dc) == 8.0
    np.testing.assert_allclose(dc.GetDataArray(), np.full((2, 2), 8.0))
    from scipy.sparse import csr_matrix

    sa = mfem.SparseMatrix(csr_matrix(np.diag([2.0, 3.0])))
    sb = mfem.SparseMatrix(csr_matrix(np.diag([4.0, 5.0])))
    from mfem.common.sparse_utils import sparsemat_to_scipycsr

    np.testing.assert_allclose(
        sparsemat_to_scipycsr(sparse_add(sa, sb)).toarray(), np.diag([16.0, 22.0])
    )
    np.testing.assert_allclose(
        sparsemat_to_scipycsr(sparse_scaled_add(sa, sb)).toarray(),
        np.diag([16.0, 21.0]),
    )
    dc.Assign(1.0)
    sparse_dense_add(sa, dc)
    np.testing.assert_allclose(dc.GetDataArray(), [[5.0, 1.0], [1.0, 7.0]])
    product = sparse_dense_product(sa, db)
    assert product.thisown
    np.testing.assert_allclose(product.GetDataArray(), [[4.0, 4.0], [6.0, 6.0]])
    for unsupported, args in [(deferred, (sa, sb)), (deferred_none, (sa, sb)),
                              (deferred_output, (sa, sb, sa))]:
        try:
            unsupported(*args)
        except TypingError:
            pass
        else:
            raise AssertionError("deferred sparse Mult compiled")
    # Product blocks are owned, but the offsets remain external: keep them live.
    offsets = mfem.intArray([0, 2])
    ba, bb = mfem.BlockMatrix(offsets), mfem.BlockMatrix(offsets)
    ba.SetBlock(0, 0, sa)
    bb.SetBlock(0, 0, sb)
    bc = block_product(ba, bb)
    assert bc.thisown
    mono = bc.CreateMonolithic()
    np.testing.assert_allclose(
        sparsemat_to_scipycsr(mono).toarray(), np.diag([8.0, 15.0])
    )
    if "-p" in sys.argv:
        from mpi4py import MPI

        mesh = mfem.Mesh(2, 2, "QUADRILATERAL")
        pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
        fec = mfem.H1_FECollection(1, mesh.Dimension())
        space = mfem.ParFiniteElementSpace(pmesh, fec)
        form = mfem.ParBilinearForm(space)
        form.AddDomainIntegrator(mfem.MassIntegrator())
        form.Assemble()
        form.Finalize()
        matrix = form.ParallelAssemble()
        owned = hypre_operations(matrix, matrix)
        assert owned.thisown
        source = mfem.Vector(matrix.Width())
        source.Assign(1.0)
        expected = mfem.Vector(matrix.Height())
        actual = mfem.Vector(matrix.Height())
        matrix.Mult(source, expected)
        owned.Mult(source, actual)
        np.testing.assert_allclose(
            actual.GetDataArray(), 12.0 * expected.GetDataArray()
        )
    print("native operators and helpers: Normal End")


if __name__ == "__main__":
    run_test()
