"""Verify the mfem::DenseTensor ordinary Numba bridge and data view."""
from numba import njit
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def data_view_round_trip(tensor):
    """Write through PyMFEM's documented [k, i, j] tensor data view."""
    data = tensor.GetDataArray()
    data[2, 1, 0] = 19.0
    return data.shape[0], data.shape[1], data.shape[2], data[2, 1, 0]


def run_test():
    tensor = mfem.DenseTensor(2, 3, 4)
    assert data_view_round_trip(tensor) == (4, 2, 3, 19.0)

    # PyMFEM maps tensor[i, j, k] to GetDataArray()[k, i, j].
    assert tensor[1, 0, 2] == 19.0


if __name__ == "__main__":
    run_test()
