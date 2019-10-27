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
    from . import _common_functions
else:
    import _common_functions

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


import mfem._ser.array
import mfem._ser.mem_manager

def InnerProduct(x, y):
    r"""InnerProduct(mfem::Vector const & x, mfem::Vector const & y) -> double"""
    return _common_functions.InnerProduct(x, y)

def RAP(*args):
    r"""
    RAP(mfem::SparseMatrix const & A, mfem::DenseMatrix & P) -> mfem::DenseMatrix
    RAP(mfem::DenseMatrix & A, mfem::SparseMatrix const & P) -> mfem::DenseMatrix
    RAP(mfem::SparseMatrix const & A, mfem::SparseMatrix const & R, mfem::SparseMatrix * ORAP=None) -> mfem::SparseMatrix
    RAP(mfem::SparseMatrix const & Rt, mfem::SparseMatrix const & A, mfem::SparseMatrix const & P) -> mfem::SparseMatrix *
    """
    return _common_functions.RAP(*args)

def Add(*args):
    r"""
    Add(mfem::DenseMatrix const & A, mfem::DenseMatrix const & B, double alpha, mfem::DenseMatrix & C)
    Add(double alpha, double const * A, double beta, double const * B, mfem::DenseMatrix & C)
    Add(double alpha, mfem::DenseMatrix const & A, double beta, mfem::DenseMatrix const & B, mfem::DenseMatrix & C)
    Add(mfem::SparseMatrix const & A, mfem::SparseMatrix const & B) -> mfem::SparseMatrix
    Add(double a, mfem::SparseMatrix const & A, double b, mfem::SparseMatrix const & B) -> mfem::SparseMatrix
    Add(mfem::Array< mfem::SparseMatrix * > & Ai) -> mfem::SparseMatrix
    Add(mfem::SparseMatrix const & A, double alpha, mfem::DenseMatrix & B)
    """
    return _common_functions.Add(*args)

def Transpose(*args):
    r"""
    Transpose(mfem::Table const & A, mfem::Table & At, int _ncols_A=-1)
    Transpose(mfem::Table const & A) -> mfem::Table
    Transpose(intArray A, mfem::Table & At, int _ncols_A=-1)
    Transpose(mfem::SparseMatrix const & A) -> mfem::SparseMatrix
    Transpose(mfem::BlockMatrix const & A) -> mfem::BlockMatrix *
    """
    return _common_functions.Transpose(*args)

def Mult(*args):
    r"""
    Mult(mfem::Table const & A, mfem::Table const & B, mfem::Table & C)
    Mult(mfem::Table const & A, mfem::Table const & B) -> mfem::Table
    Mult(mfem::DenseMatrix const & b, mfem::DenseMatrix const & c, mfem::DenseMatrix & a)
    Mult(mfem::SparseMatrix const & A, mfem::SparseMatrix const & B, mfem::SparseMatrix * OAB=None) -> mfem::SparseMatrix
    Mult(mfem::SparseMatrix const & A, mfem::DenseMatrix & B) -> mfem::DenseMatrix
    Mult(mfem::BlockMatrix const & A, mfem::BlockMatrix const & B) -> mfem::BlockMatrix *
    """
    return _common_functions.Mult(*args)


