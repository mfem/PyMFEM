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
    from . import _fespace
else:
    import _fespace

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

import mfem._ser.array
import mfem._ser.mem_manager
import mfem._ser.vector
import mfem._ser.coefficient
import mfem._ser.matrix
import mfem._ser.operators
import mfem._ser.intrules
import mfem._ser.sparsemat
import mfem._ser.densemat
import mfem._ser.eltrans
import mfem._ser.fe
import mfem._ser.geom
import mfem._ser.mesh
import mfem._ser.ncmesh
import mfem._ser.gridfunc
import mfem._ser.bilininteg
import mfem._ser.fe_coll
import mfem._ser.lininteg
import mfem._ser.linearform
import mfem._ser.element
import mfem._ser.table
import mfem._ser.hash
import mfem._ser.vertex
import mfem._ser.handle
class Ordering(object):
    r"""Proxy of C++ mfem::Ordering class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    byNODES = _fespace.Ordering_byNODES
    
    byVDIM = _fespace.Ordering_byVDIM
    

    def __init__(self):
        r"""__init__(Ordering self) -> Ordering"""
        _fespace.Ordering_swiginit(self, _fespace.new_Ordering())
    __swig_destroy__ = _fespace.delete_Ordering

# Register Ordering in _fespace:
_fespace.Ordering_swigregister(Ordering)

ElementDofOrdering_NATIVE = _fespace.ElementDofOrdering_NATIVE

ElementDofOrdering_LEXICOGRAPHIC = _fespace.ElementDofOrdering_LEXICOGRAPHIC

class FiniteElementSpace(object):
    r"""Proxy of C++ mfem::FiniteElementSpace class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(FiniteElementSpace self) -> FiniteElementSpace
        __init__(FiniteElementSpace self, FiniteElementSpace orig, Mesh mesh=None, FiniteElementCollection fec=None) -> FiniteElementSpace
        __init__(FiniteElementSpace self, Mesh mesh, FiniteElementCollection fec, int vdim=1, int ordering=byNODES) -> FiniteElementSpace
        __init__(FiniteElementSpace self, Mesh mesh, mfem::NURBSExtension * ext, FiniteElementCollection fec, int vdim=1, int ordering=byNODES) -> FiniteElementSpace
        """
        _fespace.FiniteElementSpace_swiginit(self, _fespace.new_FiniteElementSpace(*args))


        '''
            FiniteElementSpace(Mesh *m, const FiniteElementCollection *f,
                          int vdim = 1, int ordering = Ordering::byNODES);
            keep reference to mesh and fec to prevent it from
            freeed from pytho garbage collector
        '''
        if hasattr(self, 'this'):
            self.mesh = args[0]
            self.fec = args[1]





    def GetMesh(self):
        r"""GetMesh(FiniteElementSpace self) -> Mesh"""
        return _fespace.FiniteElementSpace_GetMesh(self)

    def GetNURBSext(self, *args):
        r"""
        GetNURBSext(FiniteElementSpace self) -> mfem::NURBSExtension const
        GetNURBSext(FiniteElementSpace self) -> mfem::NURBSExtension *
        """
        return _fespace.FiniteElementSpace_GetNURBSext(self, *args)

    def StealNURBSext(self):
        r"""StealNURBSext(FiniteElementSpace self) -> mfem::NURBSExtension *"""
        return _fespace.FiniteElementSpace_StealNURBSext(self)

    def Conforming(self):
        r"""Conforming(FiniteElementSpace self) -> bool"""
        return _fespace.FiniteElementSpace_Conforming(self)

    def Nonconforming(self):
        r"""Nonconforming(FiniteElementSpace self) -> bool"""
        return _fespace.FiniteElementSpace_Nonconforming(self)

    def GetConformingProlongation(self):
        r"""GetConformingProlongation(FiniteElementSpace self) -> SparseMatrix"""
        return _fespace.FiniteElementSpace_GetConformingProlongation(self)

    def GetConformingRestriction(self):
        r"""GetConformingRestriction(FiniteElementSpace self) -> SparseMatrix"""
        return _fespace.FiniteElementSpace_GetConformingRestriction(self)

    def GetProlongationMatrix(self):
        r"""GetProlongationMatrix(FiniteElementSpace self) -> Operator"""
        return _fespace.FiniteElementSpace_GetProlongationMatrix(self)

    def GetRestrictionMatrix(self):
        r"""GetRestrictionMatrix(FiniteElementSpace self) -> SparseMatrix"""
        return _fespace.FiniteElementSpace_GetRestrictionMatrix(self)

    def GetElementRestriction(self, e_ordering):
        r"""GetElementRestriction(FiniteElementSpace self, mfem::ElementDofOrdering e_ordering) -> Operator"""
        return _fespace.FiniteElementSpace_GetElementRestriction(self, e_ordering)

    def GetQuadratureInterpolator(self, *args):
        r"""
        GetQuadratureInterpolator(FiniteElementSpace self, IntegrationRule ir) -> QuadratureInterpolator
        GetQuadratureInterpolator(FiniteElementSpace self, QuadratureSpace qs) -> QuadratureInterpolator
        """
        return _fespace.FiniteElementSpace_GetQuadratureInterpolator(self, *args)

    def GetVDim(self):
        r"""GetVDim(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetVDim(self)

    def GetOrder(self, i):
        r"""GetOrder(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetOrder(self, i)

    def GetFaceOrder(self, i):
        r"""GetFaceOrder(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetFaceOrder(self, i)

    def GetNDofs(self):
        r"""GetNDofs(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNDofs(self)

    def GetVSize(self):
        r"""GetVSize(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetVSize(self)

    def GetTrueVSize(self):
        r"""GetTrueVSize(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetTrueVSize(self)

    def GetNConformingDofs(self):
        r"""GetNConformingDofs(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNConformingDofs(self)

    def GetConformingVSize(self):
        r"""GetConformingVSize(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetConformingVSize(self)

    def GetOrdering(self):
        r"""GetOrdering(FiniteElementSpace self) -> mfem::Ordering::Type"""
        return _fespace.FiniteElementSpace_GetOrdering(self)

    def FEColl(self):
        r"""FEColl(FiniteElementSpace self) -> FiniteElementCollection"""
        return _fespace.FiniteElementSpace_FEColl(self)

    def GetNVDofs(self):
        r"""GetNVDofs(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNVDofs(self)

    def GetNEDofs(self):
        r"""GetNEDofs(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNEDofs(self)

    def GetNFDofs(self):
        r"""GetNFDofs(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNFDofs(self)

    def GetNV(self):
        r"""GetNV(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNV(self)

    def GetNE(self):
        r"""GetNE(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNE(self)

    def GetNF(self):
        r"""GetNF(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNF(self)

    def GetNBE(self):
        r"""GetNBE(FiniteElementSpace self) -> int"""
        return _fespace.FiniteElementSpace_GetNBE(self)

    def GetElementType(self, i):
        r"""GetElementType(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetElementType(self, i)

    def GetElementVertices(self, i, vertices):
        r"""GetElementVertices(FiniteElementSpace self, int i, intArray vertices)"""
        return _fespace.FiniteElementSpace_GetElementVertices(self, i, vertices)

    def GetBdrElementType(self, i):
        r"""GetBdrElementType(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetBdrElementType(self, i)

    def GetElementTransformation(self, *args):
        r"""
        GetElementTransformation(FiniteElementSpace self, int i) -> ElementTransformation
        GetElementTransformation(FiniteElementSpace self, int i, IsoparametricTransformation ElTr)
        """
        return _fespace.FiniteElementSpace_GetElementTransformation(self, *args)

    def GetBdrElementTransformation(self, i):
        r"""GetBdrElementTransformation(FiniteElementSpace self, int i) -> ElementTransformation"""
        return _fespace.FiniteElementSpace_GetBdrElementTransformation(self, i)

    def GetAttribute(self, i):
        r"""GetAttribute(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetAttribute(self, i)

    def GetBdrAttribute(self, i):
        r"""GetBdrAttribute(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetBdrAttribute(self, i)

    def GetElementDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetElementDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetBdrElementDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetBdrElementDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetFaceDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetFaceDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetEdgeDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetEdgeDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetVertexDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetVertexDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetElementInteriorDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetElementInteriorDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetFaceInteriorDofs(self, i, dofs):
        r"""GetFaceInteriorDofs(FiniteElementSpace self, int i, intArray dofs)"""
        return _fespace.FiniteElementSpace_GetFaceInteriorDofs(self, i, dofs)

    def GetNumElementInteriorDofs(self, i):
        r"""GetNumElementInteriorDofs(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetNumElementInteriorDofs(self, i)

    def GetEdgeInteriorDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetEdgeInteriorDofs(self, i, vdofs)
        return vdofs.ToList()



    def DofsToVDofs(self, *args):
        r"""
        DofsToVDofs(FiniteElementSpace self, intArray dofs, int ndofs=-1)
        DofsToVDofs(FiniteElementSpace self, int vd, intArray dofs, int ndofs=-1)
        """
        return _fespace.FiniteElementSpace_DofsToVDofs(self, *args)

    def DofToVDof(self, dof, vd, ndofs=-1):
        r"""DofToVDof(FiniteElementSpace self, int dof, int vd, int ndofs=-1) -> int"""
        return _fespace.FiniteElementSpace_DofToVDof(self, dof, vd, ndofs)

    def VDofToDof(self, vdof):
        r"""VDofToDof(FiniteElementSpace self, int vdof) -> int"""
        return _fespace.FiniteElementSpace_VDofToDof(self, vdof)

    @staticmethod
    def AdjustVDofs(vdofs):
        r"""AdjustVDofs(intArray vdofs)"""
        return _fespace.FiniteElementSpace_AdjustVDofs(vdofs)

    def GetElementVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetElementVDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetBdrElementVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetBdrElementVDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetFaceVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetFaceVDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetEdgeVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetEdgeVDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetVertexVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetVertexVDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetElementInteriorVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetElementInteriorVDofs(self, i, vdofs)
        return vdofs.ToList()



    def GetEdgeInteriorVDofs(self, i):
        from  .array import intArray
        vdofs = intArray()
        _fespace.FiniteElementSpace_GetEdgeInteriorVDofs(self, i, vdofs)
        return vdofs.ToList()



    def RebuildElementToDofTable(self):
        r"""RebuildElementToDofTable(FiniteElementSpace self)"""
        return _fespace.FiniteElementSpace_RebuildElementToDofTable(self)

    def ReorderElementToDofTable(self):
        r"""ReorderElementToDofTable(FiniteElementSpace self)"""
        return _fespace.FiniteElementSpace_ReorderElementToDofTable(self)

    def BuildDofToArrays(self):
        r"""BuildDofToArrays(FiniteElementSpace self)"""
        return _fespace.FiniteElementSpace_BuildDofToArrays(self)

    def GetElementToDofTable(self):
        r"""GetElementToDofTable(FiniteElementSpace self) -> Table"""
        return _fespace.FiniteElementSpace_GetElementToDofTable(self)

    def GetBdrElementToDofTable(self):
        r"""GetBdrElementToDofTable(FiniteElementSpace self) -> Table"""
        return _fespace.FiniteElementSpace_GetBdrElementToDofTable(self)

    def GetElementForDof(self, i):
        r"""GetElementForDof(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetElementForDof(self, i)

    def GetLocalDofForDof(self, i):
        r"""GetLocalDofForDof(FiniteElementSpace self, int i) -> int"""
        return _fespace.FiniteElementSpace_GetLocalDofForDof(self, i)

    def GetFE(self, i):
        r"""GetFE(FiniteElementSpace self, int i) -> FiniteElement"""
        return _fespace.FiniteElementSpace_GetFE(self, i)

    def GetBE(self, i):
        r"""GetBE(FiniteElementSpace self, int i) -> FiniteElement"""
        return _fespace.FiniteElementSpace_GetBE(self, i)

    def GetFaceElement(self, i):
        r"""GetFaceElement(FiniteElementSpace self, int i) -> FiniteElement"""
        return _fespace.FiniteElementSpace_GetFaceElement(self, i)

    def GetEdgeElement(self, i):
        r"""GetEdgeElement(FiniteElementSpace self, int i) -> FiniteElement"""
        return _fespace.FiniteElementSpace_GetEdgeElement(self, i)

    def GetTraceElement(self, i, geom_type):
        r"""GetTraceElement(FiniteElementSpace self, int i, mfem::Geometry::Type geom_type) -> FiniteElement"""
        return _fespace.FiniteElementSpace_GetTraceElement(self, i, geom_type)

    def GetEssentialVDofs(self, bdr_attr_is_ess, ess_vdofs, component=-1):
        r"""GetEssentialVDofs(FiniteElementSpace self, intArray bdr_attr_is_ess, intArray ess_vdofs, int component=-1)"""
        return _fespace.FiniteElementSpace_GetEssentialVDofs(self, bdr_attr_is_ess, ess_vdofs, component)

    def GetEssentialTrueDofs(self, bdr_attr_is_ess, ess_tdof_list, component=-1):
        r"""GetEssentialTrueDofs(FiniteElementSpace self, intArray bdr_attr_is_ess, intArray ess_tdof_list, int component=-1)"""
        return _fespace.FiniteElementSpace_GetEssentialTrueDofs(self, bdr_attr_is_ess, ess_tdof_list, component)

    @staticmethod
    def MarkerToList(marker, list):
        r"""MarkerToList(intArray marker, intArray list)"""
        return _fespace.FiniteElementSpace_MarkerToList(marker, list)

    @staticmethod
    def ListToMarker(list, marker_size, marker, mark_val=-1):
        r"""ListToMarker(intArray list, int marker_size, intArray marker, int mark_val=-1)"""
        return _fespace.FiniteElementSpace_ListToMarker(list, marker_size, marker, mark_val)

    def ConvertToConformingVDofs(self, dofs, cdofs):
        r"""ConvertToConformingVDofs(FiniteElementSpace self, intArray dofs, intArray cdofs)"""
        return _fespace.FiniteElementSpace_ConvertToConformingVDofs(self, dofs, cdofs)

    def ConvertFromConformingVDofs(self, cdofs, dofs):
        r"""ConvertFromConformingVDofs(FiniteElementSpace self, intArray cdofs, intArray dofs)"""
        return _fespace.FiniteElementSpace_ConvertFromConformingVDofs(self, cdofs, dofs)

    def D2C_GlobalRestrictionMatrix(self, cfes):
        r"""D2C_GlobalRestrictionMatrix(FiniteElementSpace self, FiniteElementSpace cfes) -> SparseMatrix"""
        return _fespace.FiniteElementSpace_D2C_GlobalRestrictionMatrix(self, cfes)

    def D2Const_GlobalRestrictionMatrix(self, cfes):
        r"""D2Const_GlobalRestrictionMatrix(FiniteElementSpace self, FiniteElementSpace cfes) -> SparseMatrix"""
        return _fespace.FiniteElementSpace_D2Const_GlobalRestrictionMatrix(self, cfes)

    def H2L_GlobalRestrictionMatrix(self, lfes):
        r"""H2L_GlobalRestrictionMatrix(FiniteElementSpace self, FiniteElementSpace lfes) -> SparseMatrix"""
        return _fespace.FiniteElementSpace_H2L_GlobalRestrictionMatrix(self, lfes)

    def GetTransferOperator(self, coarse_fes, T):
        r"""GetTransferOperator(FiniteElementSpace self, FiniteElementSpace coarse_fes, OperatorHandle T)"""
        return _fespace.FiniteElementSpace_GetTransferOperator(self, coarse_fes, T)

    def GetTrueTransferOperator(self, coarse_fes, T):
        r"""GetTrueTransferOperator(FiniteElementSpace self, FiniteElementSpace coarse_fes, OperatorHandle T)"""
        return _fespace.FiniteElementSpace_GetTrueTransferOperator(self, coarse_fes, T)

    def Update(self, want_transform=True):
        r"""Update(FiniteElementSpace self, bool want_transform=True)"""
        return _fespace.FiniteElementSpace_Update(self, want_transform)

    def GetUpdateOperator(self, *args):
        r"""
        GetUpdateOperator(FiniteElementSpace self) -> Operator
        GetUpdateOperator(FiniteElementSpace self, OperatorHandle T)
        """
        return _fespace.FiniteElementSpace_GetUpdateOperator(self, *args)

    def SetUpdateOperatorOwner(self, own):
        r"""SetUpdateOperatorOwner(FiniteElementSpace self, bool own)"""
        return _fespace.FiniteElementSpace_SetUpdateOperatorOwner(self, own)

    def SetUpdateOperatorType(self, tid):
        r"""SetUpdateOperatorType(FiniteElementSpace self, mfem::Operator::Type tid)"""
        return _fespace.FiniteElementSpace_SetUpdateOperatorType(self, tid)

    def UpdatesFinished(self):
        r"""UpdatesFinished(FiniteElementSpace self)"""
        return _fespace.FiniteElementSpace_UpdatesFinished(self)

    def GetSequence(self):
        r"""GetSequence(FiniteElementSpace self) -> long"""
        return _fespace.FiniteElementSpace_GetSequence(self)

    def Load(self, m, input):
        r"""Load(FiniteElementSpace self, Mesh m, std::istream & input) -> FiniteElementCollection"""
        return _fespace.FiniteElementSpace_Load(self, m, input)
    __swig_destroy__ = _fespace.delete_FiniteElementSpace

    def Save(self, *args):
        r"""
        Save(FiniteElementSpace self, std::ostream & out)
        Save(FiniteElementSpace self, char const * file, int precision=8)
        Save(FiniteElementSpace self)
        """
        return _fespace.FiniteElementSpace_Save(self, *args)

# Register FiniteElementSpace in _fespace:
_fespace.FiniteElementSpace_swigregister(FiniteElementSpace)

def FiniteElementSpace_AdjustVDofs(vdofs):
    r"""FiniteElementSpace_AdjustVDofs(intArray vdofs)"""
    return _fespace.FiniteElementSpace_AdjustVDofs(vdofs)

def FiniteElementSpace_MarkerToList(marker, list):
    r"""FiniteElementSpace_MarkerToList(intArray marker, intArray list)"""
    return _fespace.FiniteElementSpace_MarkerToList(marker, list)

def FiniteElementSpace_ListToMarker(list, marker_size, marker, mark_val=-1):
    r"""FiniteElementSpace_ListToMarker(intArray list, int marker_size, intArray marker, int mark_val=-1)"""
    return _fespace.FiniteElementSpace_ListToMarker(list, marker_size, marker, mark_val)

class QuadratureSpace(object):
    r"""Proxy of C++ mfem::QuadratureSpace class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(QuadratureSpace self, Mesh mesh_, int order_) -> QuadratureSpace
        __init__(QuadratureSpace self, Mesh mesh_, std::istream & _in) -> QuadratureSpace
        """
        _fespace.QuadratureSpace_swiginit(self, _fespace.new_QuadratureSpace(*args))
    __swig_destroy__ = _fespace.delete_QuadratureSpace

    def GetSize(self):
        r"""GetSize(QuadratureSpace self) -> int"""
        return _fespace.QuadratureSpace_GetSize(self)

    def GetElementIntRule(self, idx):
        r"""GetElementIntRule(QuadratureSpace self, int idx) -> IntegrationRule"""
        return _fespace.QuadratureSpace_GetElementIntRule(self, idx)

    def Save(self, *args):
        r"""
        Save(QuadratureSpace self, std::ostream & out)
        Save(QuadratureSpace self, char const * file, int precision=8)
        Save(QuadratureSpace self)
        """
        return _fespace.QuadratureSpace_Save(self, *args)

# Register QuadratureSpace in _fespace:
_fespace.QuadratureSpace_swigregister(QuadratureSpace)

class GridTransfer(object):
    r"""Proxy of C++ mfem::GridTransfer class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _fespace.delete_GridTransfer

    def SetOperatorType(self, type):
        r"""SetOperatorType(GridTransfer self, mfem::Operator::Type type)"""
        return _fespace.GridTransfer_SetOperatorType(self, type)

    def ForwardOperator(self):
        r"""ForwardOperator(GridTransfer self) -> Operator"""
        return _fespace.GridTransfer_ForwardOperator(self)

    def BackwardOperator(self):
        r"""BackwardOperator(GridTransfer self) -> Operator"""
        return _fespace.GridTransfer_BackwardOperator(self)

    def TrueForwardOperator(self):
        r"""TrueForwardOperator(GridTransfer self) -> Operator"""
        return _fespace.GridTransfer_TrueForwardOperator(self)

    def TrueBackwardOperator(self):
        r"""TrueBackwardOperator(GridTransfer self) -> Operator"""
        return _fespace.GridTransfer_TrueBackwardOperator(self)

# Register GridTransfer in _fespace:
_fespace.GridTransfer_swigregister(GridTransfer)

class InterpolationGridTransfer(GridTransfer):
    r"""Proxy of C++ mfem::InterpolationGridTransfer class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, coarse_fes, fine_fes):
        r"""__init__(InterpolationGridTransfer self, FiniteElementSpace coarse_fes, FiniteElementSpace fine_fes) -> InterpolationGridTransfer"""
        _fespace.InterpolationGridTransfer_swiginit(self, _fespace.new_InterpolationGridTransfer(coarse_fes, fine_fes))
    __swig_destroy__ = _fespace.delete_InterpolationGridTransfer

    def SetMassIntegrator(self, mass_integ_, own_mass_integ_=True):
        r"""SetMassIntegrator(InterpolationGridTransfer self, BilinearFormIntegrator mass_integ_, bool own_mass_integ_=True)"""
        return _fespace.InterpolationGridTransfer_SetMassIntegrator(self, mass_integ_, own_mass_integ_)

    def ForwardOperator(self):
        r"""ForwardOperator(InterpolationGridTransfer self) -> Operator"""
        return _fespace.InterpolationGridTransfer_ForwardOperator(self)

    def BackwardOperator(self):
        r"""BackwardOperator(InterpolationGridTransfer self) -> Operator"""
        return _fespace.InterpolationGridTransfer_BackwardOperator(self)

# Register InterpolationGridTransfer in _fespace:
_fespace.InterpolationGridTransfer_swigregister(InterpolationGridTransfer)

class L2ProjectionGridTransfer(GridTransfer):
    r"""Proxy of C++ mfem::L2ProjectionGridTransfer class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, coarse_fes, fine_fes):
        r"""__init__(L2ProjectionGridTransfer self, FiniteElementSpace coarse_fes, FiniteElementSpace fine_fes) -> L2ProjectionGridTransfer"""
        _fespace.L2ProjectionGridTransfer_swiginit(self, _fespace.new_L2ProjectionGridTransfer(coarse_fes, fine_fes))

    def ForwardOperator(self):
        r"""ForwardOperator(L2ProjectionGridTransfer self) -> Operator"""
        return _fespace.L2ProjectionGridTransfer_ForwardOperator(self)

    def BackwardOperator(self):
        r"""BackwardOperator(L2ProjectionGridTransfer self) -> Operator"""
        return _fespace.L2ProjectionGridTransfer_BackwardOperator(self)
    __swig_destroy__ = _fespace.delete_L2ProjectionGridTransfer

# Register L2ProjectionGridTransfer in _fespace:
_fespace.L2ProjectionGridTransfer_swigregister(L2ProjectionGridTransfer)

class ElementRestriction(mfem._ser.operators.Operator):
    r"""Proxy of C++ mfem::ElementRestriction class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, arg2, arg3):
        r"""__init__(ElementRestriction self, FiniteElementSpace arg2, mfem::ElementDofOrdering arg3) -> ElementRestriction"""
        _fespace.ElementRestriction_swiginit(self, _fespace.new_ElementRestriction(arg2, arg3))

    def Mult(self, x, y):
        r"""Mult(ElementRestriction self, Vector x, Vector y)"""
        return _fespace.ElementRestriction_Mult(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(ElementRestriction self, Vector x, Vector y)"""
        return _fespace.ElementRestriction_MultTranspose(self, x, y)
    __swig_destroy__ = _fespace.delete_ElementRestriction

# Register ElementRestriction in _fespace:
_fespace.ElementRestriction_swigregister(ElementRestriction)

class QuadratureInterpolator(object):
    r"""Proxy of C++ mfem::QuadratureInterpolator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    VALUES = _fespace.QuadratureInterpolator_VALUES
    
    DERIVATIVES = _fespace.QuadratureInterpolator_DERIVATIVES
    
    DETERMINANTS = _fespace.QuadratureInterpolator_DETERMINANTS
    

    def __init__(self, *args):
        r"""
        __init__(QuadratureInterpolator self, FiniteElementSpace fes, IntegrationRule ir) -> QuadratureInterpolator
        __init__(QuadratureInterpolator self, FiniteElementSpace fes, QuadratureSpace qs) -> QuadratureInterpolator
        """
        _fespace.QuadratureInterpolator_swiginit(self, _fespace.new_QuadratureInterpolator(*args))

    def DisableTensorProducts(self, disable=True):
        r"""DisableTensorProducts(QuadratureInterpolator self, bool disable=True)"""
        return _fespace.QuadratureInterpolator_DisableTensorProducts(self, disable)

    def Mult(self, e_vec, eval_flags, q_val, q_der, q_det):
        r"""Mult(QuadratureInterpolator self, Vector e_vec, unsigned int eval_flags, Vector q_val, Vector q_der, Vector q_det)"""
        return _fespace.QuadratureInterpolator_Mult(self, e_vec, eval_flags, q_val, q_der, q_det)

    def MultTranspose(self, eval_flags, q_val, q_der, e_vec):
        r"""MultTranspose(QuadratureInterpolator self, unsigned int eval_flags, Vector q_val, Vector q_der, Vector e_vec)"""
        return _fespace.QuadratureInterpolator_MultTranspose(self, eval_flags, q_val, q_der, e_vec)
    __swig_destroy__ = _fespace.delete_QuadratureInterpolator

# Register QuadratureInterpolator in _fespace:
_fespace.QuadratureInterpolator_swigregister(QuadratureInterpolator)



