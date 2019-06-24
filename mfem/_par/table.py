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
        mname = '.'.join((pkg, '_table')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_table')
    _table = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_table', [dirname(__file__)])
        except ImportError:
            import _table
            return _table
        try:
            _mod = imp.load_module('_table', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _table = swig_import_helper()
    del swig_import_helper
else:
    import _table
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

import mfem._par.array
import mfem._par.mem_manager
class Connection(_object):
    """Proxy of C++ mfem::Connection class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Connection, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Connection, name)
    __repr__ = _swig_repr
    __swig_setmethods__["_from"] = _table.Connection__from_set
    __swig_getmethods__["_from"] = _table.Connection__from_get
    if _newclass:
        _from = _swig_property(_table.Connection__from_get, _table.Connection__from_set)
    __swig_setmethods__["to"] = _table.Connection_to_set
    __swig_getmethods__["to"] = _table.Connection_to_get
    if _newclass:
        to = _swig_property(_table.Connection_to_get, _table.Connection_to_set)

    def __init__(self, *args):
        """
        __init__(mfem::Connection self) -> Connection
        __init__(mfem::Connection self, int arg2, int to) -> Connection
        """
        this = _table.new_Connection(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __eq__(self, rhs):
        """__eq__(Connection self, Connection rhs) -> bool"""
        return _table.Connection___eq__(self, rhs)


    def __lt__(self, rhs):
        """__lt__(Connection self, Connection rhs) -> bool"""
        return _table.Connection___lt__(self, rhs)

    __swig_destroy__ = _table.delete_Connection
    __del__ = lambda self: None
Connection_swigregister = _table.Connection_swigregister
Connection_swigregister(Connection)

class Table(_object):
    """Proxy of C++ mfem::Table class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Table, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Table, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        """
        __init__(mfem::Table self) -> Table
        __init__(mfem::Table self, Table arg2) -> Table
        __init__(mfem::Table self, int dim, int connections_per_row=3) -> Table
        __init__(mfem::Table self, int dim) -> Table
        __init__(mfem::Table self, int nrows, mfem::Array< mfem::Connection > & list) -> Table
        __init__(mfem::Table self, int nrows, int * partitioning) -> Table
        """
        this = _table.new_Table(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def MakeI(self, nrows):
        """MakeI(Table self, int nrows)"""
        return _table.Table_MakeI(self, nrows)


    def AddAColumnInRow(self, r):
        """AddAColumnInRow(Table self, int r)"""
        return _table.Table_AddAColumnInRow(self, r)


    def AddColumnsInRow(self, r, ncol):
        """AddColumnsInRow(Table self, int r, int ncol)"""
        return _table.Table_AddColumnsInRow(self, r, ncol)


    def MakeJ(self):
        """MakeJ(Table self)"""
        return _table.Table_MakeJ(self)


    def AddConnection(self, r, c):
        """AddConnection(Table self, int r, int c)"""
        return _table.Table_AddConnection(self, r, c)


    def AddConnections(self, r, c, nc):
        """AddConnections(Table self, int r, int const * c, int nc)"""
        return _table.Table_AddConnections(self, r, c, nc)


    def ShiftUpI(self):
        """ShiftUpI(Table self)"""
        return _table.Table_ShiftUpI(self)


    def SetSize(self, dim, connections_per_row):
        """SetSize(Table self, int dim, int connections_per_row)"""
        return _table.Table_SetSize(self, dim, connections_per_row)


    def SetDims(self, rows, nnz):
        """SetDims(Table self, int rows, int nnz)"""
        return _table.Table_SetDims(self, rows, nnz)


    def Size(self):
        """Size(Table self) -> int"""
        return _table.Table_Size(self)


    def Size_of_connections(self):
        """Size_of_connections(Table self) -> int"""
        return _table.Table_Size_of_connections(self)


    def __call__(self, i, j):
        """__call__(Table self, int i, int j) -> int"""
        return _table.Table___call__(self, i, j)


    def RowSize(self, i):
        """RowSize(Table self, int i) -> int"""
        return _table.Table_RowSize(self, i)


    def GetRow(self, *args):
        """
        GetRow(Table self, int i, intArray row)
        GetRow(Table self, int i) -> int const
        GetRow(Table self, int i) -> int *
        """
        return _table.Table_GetRow(self, *args)


    def GetI(self, *args):
        """
        GetI(Table self) -> int
        GetI(Table self) -> int const *
        """
        return _table.Table_GetI(self, *args)


    def GetJ(self, *args):
        """
        GetJ(Table self) -> int
        GetJ(Table self) -> int const *
        """
        return _table.Table_GetJ(self, *args)


    def GetIMemory(self, *args):
        """
        GetIMemory(Table self) -> mfem::Memory< int >
        GetIMemory(Table self) -> mfem::Memory< int > const &
        """
        return _table.Table_GetIMemory(self, *args)


    def GetJMemory(self, *args):
        """
        GetJMemory(Table self) -> mfem::Memory< int >
        GetJMemory(Table self) -> mfem::Memory< int > const &
        """
        return _table.Table_GetJMemory(self, *args)


    def SortRows(self):
        """SortRows(Table self)"""
        return _table.Table_SortRows(self)


    def SetIJ(self, newI, newJ, newsize=-1):
        """
        SetIJ(Table self, int * newI, int * newJ, int newsize=-1)
        SetIJ(Table self, int * newI, int * newJ)
        """
        return _table.Table_SetIJ(self, newI, newJ, newsize)


    def Push(self, i, j):
        """Push(Table self, int i, int j) -> int"""
        return _table.Table_Push(self, i, j)


    def Finalize(self):
        """Finalize(Table self)"""
        return _table.Table_Finalize(self)


    def MakeFromList(self, nrows, list):
        """MakeFromList(Table self, int nrows, mfem::Array< mfem::Connection > const & list)"""
        return _table.Table_MakeFromList(self, nrows, list)


    def Width(self):
        """Width(Table self) -> int"""
        return _table.Table_Width(self)


    def LoseData(self):
        """LoseData(Table self)"""
        return _table.Table_LoseData(self)


    def Load(self, arg2):
        """Load(Table self, std::istream & arg2)"""
        return _table.Table_Load(self, arg2)


    def Copy(self, copy):
        """Copy(Table self, Table copy)"""
        return _table.Table_Copy(self, copy)


    def Swap(self, other):
        """Swap(Table self, Table other)"""
        return _table.Table_Swap(self, other)


    def Clear(self):
        """Clear(Table self)"""
        return _table.Table_Clear(self)


    def MemoryUsage(self):
        """MemoryUsage(Table self) -> long"""
        return _table.Table_MemoryUsage(self)

    __swig_destroy__ = _table.delete_Table
    __del__ = lambda self: None

    def GetRowList(self, i):
        """GetRowList(Table self, int i) -> PyObject *"""
        return _table.Table_GetRowList(self, i)


    def Print(self, *args):
        """
        Print(Table self, std::ostream & out, int width=4)
        Print(Table self, std::ostream & out)
        Print(Table self)
        Print(Table self, char const * file, int precision=8)
        Print(Table self, char const * file)
        """
        return _table.Table_Print(self, *args)


    def PrintMatlab(self, *args):
        """
        PrintMatlab(Table self, std::ostream & out)
        PrintMatlab(Table self, char const * file, int precision=8)
        PrintMatlab(Table self, char const * file)
        """
        return _table.Table_PrintMatlab(self, *args)


    def Save(self, *args):
        """
        Save(Table self, std::ostream & out)
        Save(Table self, char const * file, int precision=8)
        Save(Table self, char const * file)
        Save(Table self)
        """
        return _table.Table_Save(self, *args)

Table_swigregister = _table.Table_swigregister
Table_swigregister(Table)

class STable(Table):
    """Proxy of C++ mfem::STable class."""

    __swig_setmethods__ = {}
    for _s in [Table]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, STable, name, value)
    __swig_getmethods__ = {}
    for _s in [Table]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, STable, name)
    __repr__ = _swig_repr

    def __init__(self, dim, connections_per_row=3):
        """
        __init__(mfem::STable self, int dim, int connections_per_row=3) -> STable
        __init__(mfem::STable self, int dim) -> STable
        """
        this = _table.new_STable(dim, connections_per_row)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __call__(self, i, j):
        """__call__(STable self, int i, int j) -> int"""
        return _table.STable___call__(self, i, j)


    def Push(self, i, j):
        """Push(STable self, int i, int j) -> int"""
        return _table.STable_Push(self, i, j)

    __swig_destroy__ = _table.delete_STable
    __del__ = lambda self: None
STable_swigregister = _table.STable_swigregister
STable_swigregister(STable)

class DSTable(_object):
    """Proxy of C++ mfem::DSTable class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DSTable, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DSTable, name)
    __repr__ = _swig_repr

    def __init__(self, nrows):
        """__init__(mfem::DSTable self, int nrows) -> DSTable"""
        this = _table.new_DSTable(nrows)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def NumberOfRows(self):
        """NumberOfRows(DSTable self) -> int"""
        return _table.DSTable_NumberOfRows(self)


    def NumberOfEntries(self):
        """NumberOfEntries(DSTable self) -> int"""
        return _table.DSTable_NumberOfEntries(self)


    def Push(self, a, b):
        """Push(DSTable self, int a, int b) -> int"""
        return _table.DSTable_Push(self, a, b)


    def __call__(self, a, b):
        """__call__(DSTable self, int a, int b) -> int"""
        return _table.DSTable___call__(self, a, b)

    __swig_destroy__ = _table.delete_DSTable
    __del__ = lambda self: None
DSTable_swigregister = _table.DSTable_swigregister
DSTable_swigregister(DSTable)

# This file is compatible with both classic and new-style classes.

