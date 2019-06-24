# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_pumi')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_pumi')
    _pumi = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_pumi', [dirname(__file__)])
        except ImportError:
            import _pumi
            return _pumi
        try:
            _mod = imp.load_module('_pumi', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _pumi = swig_import_helper()
    del swig_import_helper
else:
    import _pumi
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

try:
    import weakref
    weakref_proxy = weakref.proxy
except __builtin__.Exception:
    weakref_proxy = lambda x: x


MFEM_VERSION = _pumi.MFEM_VERSION
MFEM_VERSION_STRING = _pumi.MFEM_VERSION_STRING
MFEM_VERSION_TYPE = _pumi.MFEM_VERSION_TYPE
MFEM_VERSION_TYPE_RELEASE = _pumi.MFEM_VERSION_TYPE_RELEASE
MFEM_VERSION_TYPE_DEVELOPMENT = _pumi.MFEM_VERSION_TYPE_DEVELOPMENT
MFEM_VERSION_MAJOR = _pumi.MFEM_VERSION_MAJOR
MFEM_VERSION_MINOR = _pumi.MFEM_VERSION_MINOR
MFEM_VERSION_PATCH = _pumi.MFEM_VERSION_PATCH
MFEM_SOURCE_DIR = _pumi.MFEM_SOURCE_DIR
MFEM_INSTALL_DIR = _pumi.MFEM_INSTALL_DIR
MFEM_TIMER_TYPE = _pumi.MFEM_TIMER_TYPE
MFEM_HYPRE_VERSION = _pumi.MFEM_HYPRE_VERSION
class intp(_object):
    """Proxy of C++ intp class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, intp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, intp, name)
    __repr__ = _swig_repr

    def __init__(self):
        """__init__(intp self) -> intp"""
        this = _pumi.new_intp()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pumi.delete_intp
    __del__ = lambda self: None

    def assign(self, value):
        """assign(intp self, int value)"""
        return _pumi.intp_assign(self, value)


    def value(self):
        """value(intp self) -> int"""
        return _pumi.intp_value(self)


    def cast(self):
        """cast(intp self) -> int *"""
        return _pumi.intp_cast(self)


    def frompointer(t):
        """frompointer(int * t) -> intp"""
        return _pumi.intp_frompointer(t)

    frompointer = staticmethod(frompointer)
intp_swigregister = _pumi.intp_swigregister
intp_swigregister(intp)

def intp_frompointer(t):
    """intp_frompointer(int * t) -> intp"""
    return _pumi.intp_frompointer(t)

class doublep(_object):
    """Proxy of C++ doublep class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, doublep, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, doublep, name)
    __repr__ = _swig_repr

    def __init__(self):
        """__init__(doublep self) -> doublep"""
        this = _pumi.new_doublep()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pumi.delete_doublep
    __del__ = lambda self: None

    def assign(self, value):
        """assign(doublep self, double value)"""
        return _pumi.doublep_assign(self, value)


    def value(self):
        """value(doublep self) -> double"""
        return _pumi.doublep_value(self)


    def cast(self):
        """cast(doublep self) -> double *"""
        return _pumi.doublep_cast(self)


    def frompointer(t):
        """frompointer(double * t) -> doublep"""
        return _pumi.doublep_frompointer(t)

    frompointer = staticmethod(frompointer)
doublep_swigregister = _pumi.doublep_swigregister
doublep_swigregister(doublep)

def doublep_frompointer(t):
    """doublep_frompointer(double * t) -> doublep"""
    return _pumi.doublep_frompointer(t)

import mfem._par.pgridfunc
import mfem._par.pfespace
import mfem._par.operators
import mfem._par.mem_manager
import mfem._par.vector
import mfem._par.array
import mfem._par.fespace
import mfem._par.coefficient
import mfem._par.matrix
import mfem._par.intrules
import mfem._par.sparsemat
import mfem._par.densemat
import mfem._par.eltrans
import mfem._par.fe
import mfem._par.geom
import mfem._par.mesh
import mfem._par.ncmesh
import mfem._par.element
import mfem._par.table
import mfem._par.hash
import mfem._par.vertex
import mfem._par.gridfunc
import mfem._par.bilininteg
import mfem._par.fe_coll
import mfem._par.lininteg
import mfem._par.linearform
import mfem._par.handle
import mfem._par.hypre
import mfem._par.pmesh
import mfem._par.pncmesh
import mfem._par.communication
import mfem._par.sets
import mfem._par.ostream_typemap

def ParMesh2ParPumiMesh(pmesh):
    """ParMesh2ParPumiMesh(ParMesh pmesh) -> mfem::ParPumiMesh *"""
    return _pumi.ParMesh2ParPumiMesh(pmesh)
# This file is compatible with both classic and new-style classes.

