"""Verify the mfem::DenseMatrix ordinary Numba bridge."""
from numba import njit
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def resize_and_measure(matrix):
    matrix.SetSize(2, 3)
    return matrix.TotalSize()


@njit
def data_view_round_trip(matrix):
    """Write through DenseMatrix's native Fortran-layout data view."""
    data = matrix.GetDataArray()
    data[1, 2] = 17.0
    return data.shape[0], data.shape[1], data[1, 2]


def run_test():
    matrix = mfem.DenseMatrix()
    assert resize_and_measure(matrix) == 6
    assert matrix.Height() == 2
    assert matrix.Width() == 3
    assert data_view_round_trip(matrix) == (2, 3, 17.0)
    assert matrix[1, 2] == 17.0


if __name__ == "__main__":
    run_test()
