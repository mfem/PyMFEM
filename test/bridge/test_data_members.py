"""Automatic public member access, including a const MFEM callback argument."""
import numpy as np
from numba import njit, types
from numba.core.errors import TypingError
import numba_swig_bridge as nsb
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def set_point(point):
    point.x = 0.25
    point.y = 0.5
    point.z = 0.75
    point.weight = 2.0
    point.index = 4
    return point.x + point.y + point.z + point.weight


@njit
def element_number(transformation):
    return transformation.ElementNo


@nsb.director(state={"unused": types.intc})
class PointCoefficient(mfem.Coefficient):
    def __init__(self):
        super().__init__()
        self.unused = 0

    @nsb.override
    def Eval(self, transformation, point):
        return 2.0 + point.x + 3.0 * point.y + transformation.ElementNo


class Reference(mfem.PyCoefficient):
    def EvalValue(self, x):
        return 2.0 + x[0] + 3.0 * x[1]


def run_test():
    point = mfem.IntegrationPoint()
    assert set_point(point) == 3.5
    assert (point.x, point.y, point.z, point.weight, point.index) == (0.25, 0.5, 0.75, 2.0, 4)
    mesh = mfem.Mesh(1, 1, "QUADRILATERAL")
    assert element_number(mesh.GetElementTransformation(0)) == 0
    fec = mfem.H1_FECollection(1, 2)
    space = mfem.FiniteElementSpace(mesh, fec)
    actual, expected = mfem.GridFunction(space), mfem.GridFunction(space)
    coefficient, reference = PointCoefficient(), Reference()
    try:
        actual.ProjectCoefficient(coefficient)
        expected.ProjectCoefficient(reference)
        np.testing.assert_allclose(actual.GetDataArray(), expected.GetDataArray(), rtol=1e-13, atol=1e-13)
    finally:
        coefficient.close()

    try:
        @nsb.director(state={"unused": types.intc})
        class InvalidCoefficient(mfem.Coefficient):
            @nsb.override
            def Eval(self, transformation, point):
                point.weight = 0.0
                return 1.0
    except TypingError:
        pass
    else:
        raise AssertionError("const IntegrationPoint callback argument became writable")


if __name__ == '__main__':
    run_test()
