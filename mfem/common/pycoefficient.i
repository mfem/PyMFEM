// Shared Python coefficient interfaces.
// User has two options:
//   1) Override Eval directly
//   2) Supply EvalValue(x) or EvalValue(x, t).
%pythoncode %{
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


class PyCoefficient(PyCoefficientBase):
    """
    Scalar coefficient with an EvalValue(x) or EvalValue(x, t) callback.
    Subclasses may instead override Eval directly.
    """
    def __init__(self):
        PyCoefficientBase.__init__(self, _pycoefficient_uses_time(self.EvalValue))

    def _EvalPy(self, x):
        return self.EvalValue(x.GetDataArray())

    def _EvalPyT(self, x, t):
        return self.EvalValue(x.GetDataArray(), t)

    def EvalValue(self, x):
        """
        Return the scalar value; subclasses may also accept time as (x, t).
        """
        raise NotImplementedError("Override Eval or EvalValue in a coefficient subclass")


class VectorPyCoefficient(VectorPyCoefficientBase):
    """
    Vector coefficient with an EvalValue(x) or EvalValue(x, t) callback.
    Subclasses may instead override Eval directly.
    """
    def __init__(self, dim):
        self.vdim = dim
        VectorPyCoefficientBase.__init__(
            self, dim, _pycoefficient_uses_time(self.EvalValue))

    def _EvalPy(self, x, V):
        V.Assign(self.EvalValue(x.GetDataArray()))

    def _EvalPyT(self, x, t, V):
        V.Assign(self.EvalValue(x.GetDataArray(), t))

    def EvalValue(self, x):
        """
        Return the vector value; subclasses may also accept time as (x, t).
        """
        raise NotImplementedError("Override Eval or EvalValue in a coefficient subclass")


class MatrixPyCoefficient(MatrixPyCoefficientBase):
    """
    Matrix coefficient with an EvalValue(x) or EvalValue(x, t) callback.
    Subclasses may instead override Eval directly.
    """
    def __init__(self, dim):
        self.vdim = dim
        MatrixPyCoefficientBase.__init__(
            self, dim, _pycoefficient_uses_time(self.EvalValue))

    def _EvalPy(self, x, K):
        K.Assign(self.EvalValue(x.GetDataArray()))

    def _EvalPyT(self, x, t, K):
        K.Assign(self.EvalValue(x.GetDataArray(), t))

    def EvalValue(self, x):
        """
        Return the matrix value; subclasses may also accept time as (x, t).
        """
        raise NotImplementedError("Override Eval or EvalValue in a coefficient subclass")
%}
