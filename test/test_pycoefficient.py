"""
Python coefficient signatures and MFEM evaluation-time propagation.
"""
import sys
import unittest
from unittest.mock import patch

import numpy as np

if "-p" in sys.argv:
    sys.argv.remove("-p")
    import mfem.par as mfem
else:
    import mfem.ser as mfem


class PythonCoefficientTest(unittest.TestCase):
    def setUp(self):
        self.mesh = mfem.Mesh.MakeCartesian2D(1, 1, mfem.Element.QUADRILATERAL)
        self.transformation = self.mesh.GetElementTransformation(0)
        self.point = mfem.IntegrationPoint()
        self.point.Set2(0.25, 0.5)
        self.transformation.SetIntPoint(self.point)

    def evaluate(self, coefficient):
        if isinstance(coefficient, mfem.MatrixCoefficient):
            result = mfem.DenseMatrix(2)
            coefficient.Eval(result, self.transformation, self.point)
            return result.GetDataArray().copy()
        if isinstance(coefficient, mfem.VectorCoefficient):
            result = mfem.Vector(2)
            coefficient.Eval(result, self.transformation, self.point)
            return result.GetDataArray().copy()
        return coefficient.Eval(self.transformation, self.point)

    def coefficients(self, callback):
        for base, shape in ((mfem.PyCoefficient, ()),
                            (mfem.VectorPyCoefficient, (2,)),
                            (mfem.MatrixPyCoefficient, (2, 2))):
            cls = type("TestCoefficient", (base,), {"EvalValue": callback(shape)})
            yield cls() if not shape else cls(2)

    def test_missing_callback_raises_on_evaluation(self):
        for base in (mfem.PyCoefficient, mfem.VectorPyCoefficient,
                     mfem.MatrixPyCoefficient):
            incomplete = type("IncompleteCoefficient", (base,), {})
            for cls in (base, incomplete):
                with self.subTest(cls=cls):
                    coefficient = cls() if base is mfem.PyCoefficient else cls(2)
                    with self.assertRaisesRegex(NotImplementedError, "Eval or EvalValue"):
                        self.evaluate(coefficient)

    def test_spatial(self):
        def callback(shape):
            def value(self, x):
                return np.full(shape, x[0] + x[1]) if shape else x[0] + x[1]
            return value

        for coefficient in self.coefficients(callback):
            for time in (0.0, 2.5, 0.0):
                coefficient.SetTime(time)
                np.testing.assert_allclose(self.evaluate(coefficient), 0.75)

    def test_time(self):
        def callback(shape):
            def value(self, x, time):
                result = x[0] + x[1] + time
                return np.full(shape, result) if shape else result
            return value

        for coefficient in self.coefficients(callback):
            # Signature inspection must happen only during construction.
            with patch("inspect.signature", side_effect=AssertionError("reinspection")):
                for time in (0.0, 2.5, -1.0, 0.0):
                    coefficient.SetTime(time)
                    np.testing.assert_allclose(self.evaluate(coefficient), 0.75 + time)

    def test_optional_time_and_inheritance(self):
        class Parent(mfem.PyCoefficient):
            def EvalValue(self, x, t=99.0):
                return x[0] + t

        class Child(Parent):
            pass

        coefficient = Child()
        self.assertAlmostEqual(self.evaluate(coefficient), 0.25)
        coefficient.SetTime(3.0)
        self.assertAlmostEqual(self.evaluate(coefficient), 3.25)

    def test_variadic_spatial(self):
        class Spatial(mfem.PyCoefficient):
            def EvalValue(self, x, *args):
                assert not args
                return x[0]

        coefficient = Spatial()
        coefficient.SetTime(3.0)
        self.assertAlmostEqual(self.evaluate(coefficient), 0.25)

    def test_invalid_signatures(self):
        def extra_argument(self, x, t, extra):
            raise AssertionError("callback must not execute at construction")

        def keyword_time(self, x, *, t):
            raise AssertionError("callback must not execute at construction")

        for base in (mfem.PyCoefficient, mfem.VectorPyCoefficient,
                     mfem.MatrixPyCoefficient):
            for callback in (extra_argument, keyword_time):
                cls = type("InvalidCoefficient", (base,), {"EvalValue": callback})
                with self.subTest(base=base, callback=callback):
                    with self.assertRaisesRegex(TypeError, "EvalValue must accept"):
                        cls() if base is mfem.PyCoefficient else cls(2)

    def test_callback_typeerror(self):
        calls = []

        def callback(shape):
            def value(self, x, t):
                calls.append(t)
                raise TypeError("callback body failed")
            return value

        for coefficient in self.coefficients(callback):
            with self.assertRaisesRegex(TypeError, "callback body failed"):
                self.evaluate(coefficient)
        self.assertEqual(calls, [0.0, 0.0, 0.0])

    def test_native_projection_and_assembly(self):
        class Scalar(mfem.PyCoefficient):
            def EvalValue(self, x, t):
                return 1.0 + t

        class Vector(mfem.VectorPyCoefficient):
            def EvalValue(self, x, t):
                return [1.0 + t, 2.0 + t]

        class Matrix(mfem.MatrixPyCoefficient):
            def EvalValue(self, x, t):
                return np.eye(2) * (1.0 + t)

        self.check_native_projection_and_assembly(Scalar(), Vector(2), Matrix(2))

    def test_eval_override_without_evalvalue(self):
        class Scalar(mfem.PyCoefficient):
            def Eval(self, transformation, point):
                return 1.0 + self.GetTime()

        class Vector(mfem.VectorPyCoefficient):
            def Eval(self, result, transformation, point):
                result.SetSize(2)
                result.Assign([1.0 + self.GetTime(), 2.0 + self.GetTime()])

        class Matrix(mfem.MatrixPyCoefficient):
            def Eval(self, result, transformation, point):
                result.SetSize(2, 2)
                result.Assign(np.eye(2) * (1.0 + self.GetTime()))

        self.check_native_projection_and_assembly(Scalar(), Vector(2), Matrix(2))

    def check_native_projection_and_assembly(self, scalar, vector, matrix):
        fec = mfem.H1_FECollection(1, 2)
        space = mfem.FiniteElementSpace(self.mesh, fec)
        vspace = mfem.FiniteElementSpace(self.mesh, fec, 2)
        field, vfield = mfem.GridFunction(space), mfem.GridFunction(vspace)
        for time in (0.0, 2.5, 0.0):
            for coefficient in (scalar, vector, matrix):
                coefficient.SetTime(time)
            field.ProjectCoefficient(scalar)
            np.testing.assert_allclose(field.GetDataArray(), 1.0 + time)
            vfield.ProjectCoefficient(vector)
            np.testing.assert_allclose(vfield.GetDataArray()[:space.GetNDofs()], 1.0 + time)
            np.testing.assert_allclose(vfield.GetDataArray()[space.GetNDofs():], 2.0 + time)
            reference = mfem.ConstantCoefficient(1.0 + time)
            forms = []
            for coefficient in (matrix, reference):
                form = mfem.BilinearForm(space)
                form.AddDomainIntegrator(mfem.DiffusionIntegrator(coefficient))
                form.Assemble()
                form.Finalize()
                forms.append(form)
            np.testing.assert_allclose(forms[0].SpMat().GetDataArray(),
                                       forms[1].SpMat().GetDataArray(), atol=1e-14)

    def test_removed_names(self):
        for name in ("PyCoefficientT", "VectorPyCoefficientT", "MatrixPyCoefficientT"):
            self.assertFalse(hasattr(mfem, name), name)


if __name__ == "__main__":
    unittest.main()
