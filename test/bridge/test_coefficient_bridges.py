"""Direct MFEM coefficient callbacks exercised by projection and assembly."""
import numpy as np
from numba import types
import numba_swig_bridge as nsb
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@nsb.director(state={"point": mfem.Vector}, fallback="silent")
class Scalar(mfem.Coefficient):
    def __init__(self):
        super().__init__()
        self.point = mfem.Vector(2)

    @nsb.override(stateaccess="point")
    def Eval(self, transformation, ip):
        transformation.Transform(ip, self.point)
        p = self.point.GetDataArray()
        return 1.0 + p[0] + 2.0 * p[1]


@nsb.director(state={"point": mfem.Vector}, fallback="silent")
class Vector(mfem.VectorCoefficient):
    def __init__(self):
        super().__init__(2)
        self.point = mfem.Vector(2)

    @nsb.override(stateaccess="point")
    def Eval(self, result, transformation, ip):
        transformation.Transform(ip, self.point)
        p = self.point.GetDataArray()
        result.SetSize(2)
        out = result.GetDataArray()
        out[0] = 1.0 + p[0]
        out[1] = 2.0 + p[1]


@nsb.director(state={"unused": types.intc}, fallback="silent")
class Matrix(mfem.MatrixCoefficient):
    def __init__(self):
        super().__init__(2)
        self.unused = 0

    @nsb.override(stateaccess="")
    def Eval(self, result, transformation, ip):
        result.SetSize(2, 2)
        out = result.GetDataArray()
        out[0, 0] = 2.0
        out[1, 0] = 0.5
        out[0, 1] = 0.25
        out[1, 1] = 3.0


class PythonScalar(mfem.PyCoefficient):
    def EvalValue(self, p):
        return 1.0 + p[0] + 2.0 * p[1]


class PythonVector(mfem.VectorPyCoefficient):
    def __init__(self):
        super().__init__(2)

    def EvalValue(self, p):
        return np.array([1.0 + p[0], 2.0 + p[1]])


def run_test():
    mesh = mfem.Mesh(2, 2, "QUADRILATERAL")
    fec = mfem.H1_FECollection(1, 2)
    space = mfem.FiniteElementSpace(mesh, fec)
    vspace = mfem.FiniteElementSpace(mesh, fec, 2)
    scalar, vector, matrix = Scalar(), Vector(), Matrix()
    try:
        for fs, native, reference in [(space, scalar, PythonScalar()),
                                      (vspace, vector, PythonVector())]:
            actual, expected = mfem.GridFunction(fs), mfem.GridFunction(fs)
            actual.ProjectCoefficient(native)
            expected.ProjectCoefficient(reference)
            np.testing.assert_allclose(actual.GetDataArray(), expected.GetDataArray(),
                                       rtol=1e-13, atol=1e-13)
        values = mfem.DenseMatrix(2)
        values[0, 0], values[1, 0] = 2.0, 0.5
        values[0, 1], values[1, 1] = 0.25, 3.0
        reference = mfem.MatrixConstantCoefficient(values)
        forms = []
        for coefficient in (matrix, reference):
            form = mfem.BilinearForm(space)
            form.AddDomainIntegrator(mfem.DiffusionIntegrator(coefficient))
            form.Assemble()
            form.Finalize()
            forms.append(form)
        np.testing.assert_allclose(forms[0].SpMat().GetDataArray(),
                                   forms[1].SpMat().GetDataArray(), rtol=1e-13, atol=1e-13)
    finally:
        scalar.close()
        vector.close()
        matrix.close()


if __name__ == "__main__":
    run_test()
