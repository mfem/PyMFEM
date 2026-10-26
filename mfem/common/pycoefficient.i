// Shared Python callback interfaces for coefficient class
// User should provilde either EvalValue(x) or EvalValue(x, t)
// Do not provide both. Number of arguments is used to check if
// it is time-dependent coefficient or not.
%pythoncode %{
from abc import ABC as _ABC, abstractmethod as _abstractmethod


def _pycoefficient_uses_time(callback):
    """
    Select and validate the callback signature once, without evaluating it.
    """
    import inspect

    signature = inspect.signature(callback)
    positional = [parameter for parameter in signature.parameters.values()
                  if parameter.kind in (inspect.Parameter.POSITIONAL_ONLY,
                                        inspect.Parameter.POSITIONAL_OR_KEYWORD)]
    uses_time = len(positional) >= 2
    try:
        signature.bind(None, None) if uses_time else signature.bind(None)
    except TypeError as error:
        raise TypeError("EvalValue must accept (x) or (x, t)") from error
    return int(uses_time)


class PyCoefficient(PyCoefficientBase, _ABC):
    """
    Scalar coefficient with an EvalValue(x) or EvalValue(x, t) callback.

    """
    def __init__(self):
        PyCoefficientBase.__init__(self, _pycoefficient_uses_time(self.EvalValue))

    def _EvalPy(self, x):
        return self.EvalValue(x.GetDataArray())

    def _EvalPyT(self, x, t):
        return self.EvalValue(x.GetDataArray(), t)

    @_abstractmethod
    def EvalValue(self, x):
        """
        Return the scalar value; subclasses may also accept time as (x, t).
        """
        raise NotImplementedError


class VectorPyCoefficient(VectorPyCoefficientBase, _ABC):
    """
    Vector coefficient with an EvalValue(x) or EvalValue(x, t) callback.
    """
    def __init__(self, dim):
        self.vdim = dim
        VectorPyCoefficientBase.__init__(
            self, dim, _pycoefficient_uses_time(self.EvalValue))

    def _EvalPy(self, x, V):
        V.Assign(self.EvalValue(x.GetDataArray()))

    def _EvalPyT(self, x, t, V):
        V.Assign(self.EvalValue(x.GetDataArray(), t))

    @_abstractmethod
    def EvalValue(self, x):
        """
        Return the vector value; subclasses may also accept time as (x, t).
        """
        raise NotImplementedError


class MatrixPyCoefficient(MatrixPyCoefficientBase, _ABC):
    """
    Matrix coefficient with an EvalValue(x) or EvalValue(x, t) callback.
    """
    def __init__(self, dim):
        self.vdim = dim
        MatrixPyCoefficientBase.__init__(
            self, dim, _pycoefficient_uses_time(self.EvalValue))

    def _EvalPy(self, x, K):
        K.Assign(self.EvalValue(x.GetDataArray()))

    def _EvalPyT(self, x, t, K):
        K.Assign(self.EvalValue(x.GetDataArray(), t))

    @_abstractmethod
    def EvalValue(self, x):
        """
        Return the matrix value; subclasses may also accept time as (x, t).
        """
        raise NotImplementedError
%}
