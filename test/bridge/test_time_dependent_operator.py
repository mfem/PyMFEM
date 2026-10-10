"""Verify native C++ dispatch to a void-return Numba director override."""
from numba import njit, types
import numba_swig_bridge as nsb
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@nsb.director(state={"unused": types.intc}, fallback="silent")
class DoubleOperator(mfem.TimeDependentOperator):
    def __init__(self, size):
        super().__init__(size)
        self.unused = 0

    @nsb.override
    def Mult(self, source, destination):
        source_data = source.GetDataArray()
        destination_data = destination.GetDataArray()
        for index in range(source_data.shape[0]):
            destination_data[index] = 2.0 * source_data[index]


@njit
def apply(operator, source, destination):
    operator.Mult(source, destination)


def run_test():
    operator = DoubleOperator(3)
    source = mfem.Vector(3)
    destination = mfem.Vector(3)
    source.Assign(2.0)
    try:
        apply(operator, source, destination)
        assert tuple(destination.GetDataArray()) == (4.0, 4.0, 4.0)
    finally:
        operator.close()


if __name__ == "__main__":
    run_test()
