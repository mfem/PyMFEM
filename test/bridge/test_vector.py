"""Verify the mfem::Vector ordinary Numba bridge."""
from numba import njit
import mfem.ser as mfem

# Install the registrations emitted by a bridge-enabled PyMFEM build.
registration = mfem.get_bridge_registration()


@njit
def repeated_norm(vector, count):
    total = 0.0
    for _ in range(count):
        total += vector.Norml2()
    return total


def repeated_norm_python(vector, count):
    total = 0.0
    for _ in range(count):
        total += vector.Norml2()
    return total


@njit
def data_view_round_trip(vector):
    data = vector.GetDataArray()
    data[0] = 7.0
    return data[0] + data.shape[0]


def run_test():
    count = 1_000
    vector = mfem.Vector(4)
    vector.Assign(2.0)
    expected = 4.0 * count

    # This call both warms up compilation and proves that Vector is accepted by
    # Numba through the registered ordinary bridge type.
    assert repeated_norm(vector, 1) == 4.0
    assert repeated_norm(vector, count) == expected
    assert repeated_norm_python(vector, count) == expected
    assert data_view_round_trip(vector) == 11.0
    assert vector[0] == 7.0


if __name__ == "__main__":
    run_test()
