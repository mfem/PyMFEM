# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _operators
else:
    import _operators

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


import weakref

import mfem._ser.mem_manager
import mfem._ser.vector
import mfem._ser.array
class Operator(object):
    r"""Proxy of C++ mfem::Operator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(Operator self, int s=0) -> Operator
        __init__(Operator self, int h, int w) -> Operator
        """
        if self.__class__ == Operator:
            _self = None
        else:
            _self = self
        _operators.Operator_swiginit(self, _operators.new_Operator(_self, *args))

    def Height(self):
        r"""Height(Operator self) -> int"""
        return _operators.Operator_Height(self)

    def NumRows(self):
        r"""NumRows(Operator self) -> int"""
        return _operators.Operator_NumRows(self)

    def Width(self):
        r"""Width(Operator self) -> int"""
        return _operators.Operator_Width(self)

    def NumCols(self):
        r"""NumCols(Operator self) -> int"""
        return _operators.Operator_NumCols(self)

    def GetMemoryClass(self):
        r"""GetMemoryClass(Operator self) -> mfem::MemoryClass"""
        return _operators.Operator_GetMemoryClass(self)

    def Mult(self, x, y):
        r"""Mult(Operator self, Vector x, Vector y)"""
        return _operators.Operator_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(Operator self, Vector x, Vector y)"""
        return _operators.Operator_MultTranspose(self, x, y)

    def GetGradient(self, x):
        r"""GetGradient(Operator self, Vector x) -> Operator"""
        return _operators.Operator_GetGradient(self, x)

    def GetProlongation(self):
        r"""GetProlongation(Operator self) -> Operator"""
        return _operators.Operator_GetProlongation(self)

    def GetRestriction(self):
        r"""GetRestriction(Operator self) -> Operator"""
        return _operators.Operator_GetRestriction(self)

    def FormLinearSystem(self, ess_tdof_list, x, b, A, X, B, copy_interior=0):
        r"""FormLinearSystem(Operator self, intArray ess_tdof_list, Vector x, Vector b, mfem::Operator *& A, Vector X, Vector B, int copy_interior=0)"""
        return _operators.Operator_FormLinearSystem(self, ess_tdof_list, x, b, A, X, B, copy_interior)

    def RecoverFEMSolution(self, X, b, x):
        r"""RecoverFEMSolution(Operator self, Vector X, Vector b, Vector x)"""
        return _operators.Operator_RecoverFEMSolution(self, X, b, x)
    __swig_destroy__ = _operators.delete_Operator
    ANY_TYPE = _operators.Operator_ANY_TYPE
    
    MFEM_SPARSEMAT = _operators.Operator_MFEM_SPARSEMAT
    
    Hypre_ParCSR = _operators.Operator_Hypre_ParCSR
    
    PETSC_MATAIJ = _operators.Operator_PETSC_MATAIJ
    
    PETSC_MATIS = _operators.Operator_PETSC_MATIS
    
    PETSC_MATSHELL = _operators.Operator_PETSC_MATSHELL
    
    PETSC_MATNEST = _operators.Operator_PETSC_MATNEST
    
    PETSC_MATHYPRE = _operators.Operator_PETSC_MATHYPRE
    
    PETSC_MATGENERIC = _operators.Operator_PETSC_MATGENERIC
    

    def GetType(self):
        r"""GetType(Operator self) -> mfem::Operator::Type"""
        return _operators.Operator_GetType(self)

    def PrintMatlab(self, *args):
        r"""
        PrintMatlab(Operator self, std::ostream & out, int n=0, int m=0)
        PrintMatlab(Operator self, char const * file, int precision=8)
        """
        return _operators.Operator_PrintMatlab(self, *args)
    def __disown__(self):
        self.this.disown()
        _operators.disown_Operator(self)
        return weakref.proxy(self)

# Register Operator in _operators:
_operators.Operator_swigregister(Operator)

class TimeDependentOperator(Operator):
    r"""Proxy of C++ mfem::TimeDependentOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    EXPLICIT = _operators.TimeDependentOperator_EXPLICIT
    
    IMPLICIT = _operators.TimeDependentOperator_IMPLICIT
    
    HOMOGENEOUS = _operators.TimeDependentOperator_HOMOGENEOUS
    

    def __init__(self, *args):
        r"""
        __init__(TimeDependentOperator self, int n=0, double t_=0.0, mfem::TimeDependentOperator::Type type_=EXPLICIT) -> TimeDependentOperator
        __init__(TimeDependentOperator self, int h, int w, double t_=0.0, mfem::TimeDependentOperator::Type type_=EXPLICIT) -> TimeDependentOperator
        """
        if self.__class__ == TimeDependentOperator:
            _self = None
        else:
            _self = self
        _operators.TimeDependentOperator_swiginit(self, _operators.new_TimeDependentOperator(_self, *args))

    def GetTime(self):
        r"""GetTime(TimeDependentOperator self) -> double"""
        return _operators.TimeDependentOperator_GetTime(self)

    def SetTime(self, _t):
        r"""SetTime(TimeDependentOperator self, double const _t)"""
        return _operators.TimeDependentOperator_SetTime(self, _t)

    def isExplicit(self):
        r"""isExplicit(TimeDependentOperator self) -> bool"""
        return _operators.TimeDependentOperator_isExplicit(self)

    def isImplicit(self):
        r"""isImplicit(TimeDependentOperator self) -> bool"""
        return _operators.TimeDependentOperator_isImplicit(self)

    def isHomogeneous(self):
        r"""isHomogeneous(TimeDependentOperator self) -> bool"""
        return _operators.TimeDependentOperator_isHomogeneous(self)

    def ExplicitMult(self, x, y):
        r"""ExplicitMult(TimeDependentOperator self, Vector x, Vector y)"""
        return _operators.TimeDependentOperator_ExplicitMult(self, x, y)

    def ImplicitMult(self, x, k, y):
        r"""ImplicitMult(TimeDependentOperator self, Vector x, Vector k, Vector y)"""
        return _operators.TimeDependentOperator_ImplicitMult(self, x, k, y)

    def Mult(self, x, y):
        r"""Mult(TimeDependentOperator self, Vector x, Vector y)"""
        return _operators.TimeDependentOperator_Mult(self, x, y)

    def ImplicitSolve(self, dt, x, k):
        r"""ImplicitSolve(TimeDependentOperator self, double const dt, Vector x, Vector k)"""
        return _operators.TimeDependentOperator_ImplicitSolve(self, dt, x, k)

    def GetImplicitGradient(self, x, k, shift):
        r"""GetImplicitGradient(TimeDependentOperator self, Vector x, Vector k, double shift) -> Operator"""
        return _operators.TimeDependentOperator_GetImplicitGradient(self, x, k, shift)

    def GetExplicitGradient(self, x):
        r"""GetExplicitGradient(TimeDependentOperator self, Vector x) -> Operator"""
        return _operators.TimeDependentOperator_GetExplicitGradient(self, x)
    __swig_destroy__ = _operators.delete_TimeDependentOperator
    def __disown__(self):
        self.this.disown()
        _operators.disown_TimeDependentOperator(self)
        return weakref.proxy(self)

# Register TimeDependentOperator in _operators:
_operators.TimeDependentOperator_swigregister(TimeDependentOperator)

class Solver(Operator):
    r"""Proxy of C++ mfem::Solver class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    iterative_mode = property(_operators.Solver_iterative_mode_get, _operators.Solver_iterative_mode_set, doc=r"""iterative_mode : bool""")

    def __init__(self, *args):
        r"""
        __init__(Solver self, int s=0, bool iter_mode=False) -> Solver
        __init__(Solver self, int h, int w, bool iter_mode=False) -> Solver
        """
        if self.__class__ == Solver:
            _self = None
        else:
            _self = self
        _operators.Solver_swiginit(self, _operators.new_Solver(_self, *args))

    def SetOperator(self, op):
        r"""SetOperator(Solver self, Operator op)"""
        return _operators.Solver_SetOperator(self, op)
    __swig_destroy__ = _operators.delete_Solver
    def __disown__(self):
        self.this.disown()
        _operators.disown_Solver(self)
        return weakref.proxy(self)

# Register Solver in _operators:
_operators.Solver_swigregister(Solver)

class IdentityOperator(Operator):
    r"""Proxy of C++ mfem::IdentityOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, n):
        r"""__init__(IdentityOperator self, int n) -> IdentityOperator"""
        _operators.IdentityOperator_swiginit(self, _operators.new_IdentityOperator(n))

    def Mult(self, x, y):
        r"""Mult(IdentityOperator self, Vector x, Vector y)"""
        return _operators.IdentityOperator_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(IdentityOperator self, Vector x, Vector y)"""
        return _operators.IdentityOperator_MultTranspose(self, x, y)
    __swig_destroy__ = _operators.delete_IdentityOperator

# Register IdentityOperator in _operators:
_operators.IdentityOperator_swigregister(IdentityOperator)

class TransposeOperator(Operator):
    r"""Proxy of C++ mfem::TransposeOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(TransposeOperator self, Operator a) -> TransposeOperator
        __init__(TransposeOperator self, Operator a) -> TransposeOperator
        """
        _operators.TransposeOperator_swiginit(self, _operators.new_TransposeOperator(*args))

    def Mult(self, x, y):
        r"""Mult(TransposeOperator self, Vector x, Vector y)"""
        return _operators.TransposeOperator_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(TransposeOperator self, Vector x, Vector y)"""
        return _operators.TransposeOperator_MultTranspose(self, x, y)
    __swig_destroy__ = _operators.delete_TransposeOperator

# Register TransposeOperator in _operators:
_operators.TransposeOperator_swigregister(TransposeOperator)

class ProductOperator(Operator):
    r"""Proxy of C++ mfem::ProductOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, A, B, ownA, ownB):
        r"""__init__(ProductOperator self, Operator A, Operator B, bool ownA, bool ownB) -> ProductOperator"""
        _operators.ProductOperator_swiginit(self, _operators.new_ProductOperator(A, B, ownA, ownB))

    def Mult(self, x, y):
        r"""Mult(ProductOperator self, Vector x, Vector y)"""
        return _operators.ProductOperator_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(ProductOperator self, Vector x, Vector y)"""
        return _operators.ProductOperator_MultTranspose(self, x, y)
    __swig_destroy__ = _operators.delete_ProductOperator

# Register ProductOperator in _operators:
_operators.ProductOperator_swigregister(ProductOperator)

class RAPOperator(Operator):
    r"""Proxy of C++ mfem::RAPOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, Rt_, A_, P_):
        r"""__init__(RAPOperator self, Operator Rt_, Operator A_, Operator P_) -> RAPOperator"""
        _operators.RAPOperator_swiginit(self, _operators.new_RAPOperator(Rt_, A_, P_))

    def GetMemoryClass(self):
        r"""GetMemoryClass(RAPOperator self) -> mfem::MemoryClass"""
        return _operators.RAPOperator_GetMemoryClass(self)

    def Mult(self, x, y):
        r"""Mult(RAPOperator self, Vector x, Vector y)"""
        return _operators.RAPOperator_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(RAPOperator self, Vector x, Vector y)"""
        return _operators.RAPOperator_MultTranspose(self, x, y)
    __swig_destroy__ = _operators.delete_RAPOperator

# Register RAPOperator in _operators:
_operators.RAPOperator_swigregister(RAPOperator)

class TripleProductOperator(Operator):
    r"""Proxy of C++ mfem::TripleProductOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, A, B, C, ownA, ownB, ownC):
        r"""__init__(TripleProductOperator self, Operator A, Operator B, Operator C, bool ownA, bool ownB, bool ownC) -> TripleProductOperator"""
        _operators.TripleProductOperator_swiginit(self, _operators.new_TripleProductOperator(A, B, C, ownA, ownB, ownC))

    def GetMemoryClass(self):
        r"""GetMemoryClass(TripleProductOperator self) -> mfem::MemoryClass"""
        return _operators.TripleProductOperator_GetMemoryClass(self)

    def Mult(self, x, y):
        r"""Mult(TripleProductOperator self, Vector x, Vector y)"""
        return _operators.TripleProductOperator_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(TripleProductOperator self, Vector x, Vector y)"""
        return _operators.TripleProductOperator_MultTranspose(self, x, y)
    __swig_destroy__ = _operators.delete_TripleProductOperator

# Register TripleProductOperator in _operators:
_operators.TripleProductOperator_swigregister(TripleProductOperator)

class ConstrainedOperator(Operator):
    r"""Proxy of C++ mfem::ConstrainedOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, A, list, own_A=False):
        r"""__init__(ConstrainedOperator self, Operator A, intArray list, bool own_A=False) -> ConstrainedOperator"""
        _operators.ConstrainedOperator_swiginit(self, _operators.new_ConstrainedOperator(A, list, own_A))

    def GetMemoryClass(self):
        r"""GetMemoryClass(ConstrainedOperator self) -> mfem::MemoryClass"""
        return _operators.ConstrainedOperator_GetMemoryClass(self)

    def EliminateRHS(self, x, b):
        r"""EliminateRHS(ConstrainedOperator self, Vector x, Vector b)"""
        return _operators.ConstrainedOperator_EliminateRHS(self, x, b)

    def Mult(self, x, y):
        r"""Mult(ConstrainedOperator self, Vector x, Vector y)"""
        return _operators.ConstrainedOperator_Mult(self, x, y)
    __swig_destroy__ = _operators.delete_ConstrainedOperator

# Register ConstrainedOperator in _operators:
_operators.ConstrainedOperator_swigregister(ConstrainedOperator)

class PyOperatorBase(Operator):
    r"""Proxy of C++ mfem::PyOperatorBase class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(PyOperatorBase self, int s=0) -> PyOperatorBase
        __init__(PyOperatorBase self, int h, int w) -> PyOperatorBase
        """
        if self.__class__ == PyOperatorBase:
            _self = None
        else:
            _self = self
        _operators.PyOperatorBase_swiginit(self, _operators.new_PyOperatorBase(_self, *args))

    def Mult(self, x, y):
        r"""Mult(PyOperatorBase self, Vector x, Vector y)"""
        return _operators.PyOperatorBase_Mult(self, x, y)

    def _EvalMult(self, arg0):
        r"""_EvalMult(PyOperatorBase self, Vector arg0) -> Vector"""
        return _operators.PyOperatorBase__EvalMult(self, arg0)
    __swig_destroy__ = _operators.delete_PyOperatorBase
    def __disown__(self):
        self.this.disown()
        _operators.disown_PyOperatorBase(self)
        return weakref.proxy(self)

# Register PyOperatorBase in _operators:
_operators.PyOperatorBase_swigregister(PyOperatorBase)

class PyTimeDependentOperatorBase(TimeDependentOperator):
    r"""Proxy of C++ mfem::PyTimeDependentOperatorBase class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(PyTimeDependentOperatorBase self, int n=0, double _t=0.0) -> PyTimeDependentOperatorBase
        __init__(PyTimeDependentOperatorBase self, int h, int w, double _t=0.0) -> PyTimeDependentOperatorBase
        """
        if self.__class__ == PyTimeDependentOperatorBase:
            _self = None
        else:
            _self = self
        _operators.PyTimeDependentOperatorBase_swiginit(self, _operators.new_PyTimeDependentOperatorBase(_self, *args))

    def Mult(self, x, y):
        r"""Mult(PyTimeDependentOperatorBase self, Vector x, Vector y)"""
        return _operators.PyTimeDependentOperatorBase_Mult(self, x, y)

    def _EvalMult(self, arg0):
        r"""_EvalMult(PyTimeDependentOperatorBase self, Vector arg0) -> Vector"""
        return _operators.PyTimeDependentOperatorBase__EvalMult(self, arg0)
    __swig_destroy__ = _operators.delete_PyTimeDependentOperatorBase
    def __disown__(self):
        self.this.disown()
        _operators.disown_PyTimeDependentOperatorBase(self)
        return weakref.proxy(self)

# Register PyTimeDependentOperatorBase in _operators:
_operators.PyTimeDependentOperatorBase_swigregister(PyTimeDependentOperatorBase)


class PyOperator(PyOperatorBase):
   def __init__(self, *args):
       PyOperatorBase.__init__(self, *args)
   def _EvalMult(self, x):
       return self.EvalMult(x.GetDataArray())
   def EvalMult(self, x):
       raise NotImplementedError('you must specify this method')

class PyTimeDependentOperator(PyTimeDependentOperatorBase):
   def __init__(self, *args):  
       PyTimeDependentOperatorBase.__init__(self, *args)
   def _EvalMult(self, x):
       return self.EvalMult(x.GetDataArray())
   def EvalMult(self, x):
       raise NotImplementedError('you must specify this method')




