"""Verify the mfem::Array<int> ordinary Numba bridge."""
from numba import njit
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def array_summary(values):
    return values.Size(), values.Sum(), values.Min(), values.Max()


def run_test():
    values = mfem.intArray([4, -2, 7])
    assert array_summary(values) == (3, 9, -2, 7)


if __name__ == "__main__":
    run_test()
