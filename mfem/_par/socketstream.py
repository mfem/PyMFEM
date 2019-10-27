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
    from . import _socketstream
else:
    import _socketstream

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

import mfem._par.mesh
import mfem._par.matrix
import mfem._par.vector
import mfem._par.array
import mfem._par.mem_manager
import mfem._par.operators
import mfem._par.ncmesh
import mfem._par.element
import mfem._par.densemat
import mfem._par.geom
import mfem._par.intrules
import mfem._par.table
import mfem._par.hash
import mfem._par.vertex
import mfem._par.gridfunc
import mfem._par.coefficient
import mfem._par.sparsemat
import mfem._par.eltrans
import mfem._par.fe
import mfem._par.fespace
import mfem._par.fe_coll
import mfem._par.lininteg
import mfem._par.handle
import mfem._par.hypre
import mfem._par.bilininteg
import mfem._par.linearform
class socketbuf(object):
    r"""Proxy of C++ mfem::socketbuf class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(socketbuf self) -> socketbuf
        __init__(socketbuf self, int sd) -> socketbuf
        __init__(socketbuf self, char const [] hostname, int port) -> socketbuf
        """
        _socketstream.socketbuf_swiginit(self, _socketstream.new_socketbuf(*args))

    def attach(self, sd):
        r"""attach(socketbuf self, int sd) -> int"""
        return _socketstream.socketbuf_attach(self, sd)

    def detach(self):
        r"""detach(socketbuf self) -> int"""
        return _socketstream.socketbuf_detach(self)

    def open(self, hostname, port):
        r"""open(socketbuf self, char const [] hostname, int port) -> int"""
        return _socketstream.socketbuf_open(self, hostname, port)

    def close(self):
        r"""close(socketbuf self) -> int"""
        return _socketstream.socketbuf_close(self)

    def getsocketdescriptor(self):
        r"""getsocketdescriptor(socketbuf self) -> int"""
        return _socketstream.socketbuf_getsocketdescriptor(self)

    def is_open(self):
        r"""is_open(socketbuf self) -> bool"""
        return _socketstream.socketbuf_is_open(self)
    __swig_destroy__ = _socketstream.delete_socketbuf

# Register socketbuf in _socketstream:
_socketstream.socketbuf_swigregister(socketbuf)

class socketstream(object):
    r"""Proxy of C++ mfem::socketstream class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    secure_default = _socketstream.socketstream_secure_default
    

    def __init__(self, *args):
        r"""
        __init__(socketstream self, bool secure=secure_default) -> socketstream
        __init__(socketstream self, socketbuf buf) -> socketstream
        __init__(socketstream self, int s, bool secure=secure_default) -> socketstream
        __init__(socketstream self, char const [] hostname, int port, bool secure=secure_default) -> socketstream
        """
        _socketstream.socketstream_swiginit(self, _socketstream.new_socketstream(*args))

    def rdbuf(self):
        r"""rdbuf(socketstream self) -> socketbuf"""
        return _socketstream.socketstream_rdbuf(self)

    def open(self, hostname, port):
        r"""open(socketstream self, char const [] hostname, int port) -> int"""
        return _socketstream.socketstream_open(self, hostname, port)

    def close(self):
        r"""close(socketstream self) -> int"""
        return _socketstream.socketstream_close(self)

    def is_open(self):
        r"""is_open(socketstream self) -> bool"""
        return _socketstream.socketstream_is_open(self)
    __swig_destroy__ = _socketstream.delete_socketstream

    def precision(self, *args):
        r"""
        precision(socketstream self, int const p) -> int
        precision(socketstream self) -> int
        """
        return _socketstream.socketstream_precision(self, *args)

    def send_solution(self, mesh, gf):
        r"""send_solution(socketstream self, Mesh mesh, GridFunction gf)"""
        return _socketstream.socketstream_send_solution(self, mesh, gf)

    def send_text(self, ostr):
        r"""send_text(socketstream self, char const [] ostr)"""
        return _socketstream.socketstream_send_text(self, ostr)

    def flush(self):
        r"""flush(socketstream self)"""
        return _socketstream.socketstream_flush(self)

    def __lshift__(self, *args):
        r"""
        __lshift__(socketstream self, char const [] ostr) -> socketstream
        __lshift__(socketstream self, int const x) -> socketstream
        __lshift__(socketstream self, Mesh mesh) -> socketstream
        __lshift__(socketstream self, GridFunction gf) -> socketstream
        """
        return _socketstream.socketstream___lshift__(self, *args)

    def endline(self):
        r"""endline(socketstream self) -> socketstream"""
        return _socketstream.socketstream_endline(self)

# Register socketstream in _socketstream:
_socketstream.socketstream_swigregister(socketstream)

class socketserver(object):
    r"""Proxy of C++ mfem::socketserver class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, port, backlog=4):
        r"""__init__(socketserver self, int port, int backlog=4) -> socketserver"""
        _socketstream.socketserver_swiginit(self, _socketstream.new_socketserver(port, backlog))

    def good(self):
        r"""good(socketserver self) -> bool"""
        return _socketstream.socketserver_good(self)

    def close(self):
        r"""close(socketserver self) -> int"""
        return _socketstream.socketserver_close(self)

    def accept(self, *args):
        r"""
        accept(socketserver self) -> int
        accept(socketserver self, socketstream sockstr) -> int
        """
        return _socketstream.socketserver_accept(self, *args)
    __swig_destroy__ = _socketstream.delete_socketserver

# Register socketserver in _socketstream:
_socketstream.socketserver_swigregister(socketserver)



