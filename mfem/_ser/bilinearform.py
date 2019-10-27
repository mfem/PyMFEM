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
    from . import _bilinearform
else:
    import _bilinearform

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

class intp(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self):
        _bilinearform.intp_swiginit(self, _bilinearform.new_intp())
    __swig_destroy__ = _bilinearform.delete_intp

    def assign(self, value):
        return _bilinearform.intp_assign(self, value)

    def value(self):
        return _bilinearform.intp_value(self)

    def cast(self):
        return _bilinearform.intp_cast(self)

    @staticmethod
    def frompointer(t):
        return _bilinearform.intp_frompointer(t)

# Register intp in _bilinearform:
_bilinearform.intp_swigregister(intp)

def intp_frompointer(t):
    return _bilinearform.intp_frompointer(t)

class doublep(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self):
        _bilinearform.doublep_swiginit(self, _bilinearform.new_doublep())
    __swig_destroy__ = _bilinearform.delete_doublep

    def assign(self, value):
        return _bilinearform.doublep_assign(self, value)

    def value(self):
        return _bilinearform.doublep_value(self)

    def cast(self):
        return _bilinearform.doublep_cast(self)

    @staticmethod
    def frompointer(t):
        return _bilinearform.doublep_frompointer(t)

# Register doublep in _bilinearform:
_bilinearform.doublep_swigregister(doublep)

def doublep_frompointer(t):
    return _bilinearform.doublep_frompointer(t)

import mfem._ser.mem_manager
import mfem._ser.array
import mfem._ser.fespace
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
AssemblyLevel_FULL = _bilinearform.AssemblyLevel_FULL

AssemblyLevel_ELEMENT = _bilinearform.AssemblyLevel_ELEMENT

AssemblyLevel_PARTIAL = _bilinearform.AssemblyLevel_PARTIAL

AssemblyLevel_NONE = _bilinearform.AssemblyLevel_NONE

class BilinearForm(mfem._ser.matrix.Matrix):
    r"""Proxy of C++ mfem::BilinearForm class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(BilinearForm self) -> BilinearForm
        __init__(BilinearForm self, FiniteElementSpace f) -> BilinearForm
        __init__(BilinearForm self, FiniteElementSpace f, BilinearForm bf, int ps=0) -> BilinearForm
        """
        if self.__class__ == BilinearForm:
            _self = None
        else:
            _self = self
        _bilinearform.BilinearForm_swiginit(self, _bilinearform.new_BilinearForm(_self, *args))

    def Size(self):
        r"""Size(BilinearForm self) -> int"""
        return _bilinearform.BilinearForm_Size(self)

    def SetAssemblyLevel(self, assembly_level):
        r"""SetAssemblyLevel(BilinearForm self, mfem::AssemblyLevel assembly_level)"""
        return _bilinearform.BilinearForm_SetAssemblyLevel(self, assembly_level)

    def EnableStaticCondensation(self):
        r"""EnableStaticCondensation(BilinearForm self)"""
        return _bilinearform.BilinearForm_EnableStaticCondensation(self)

    def StaticCondensationIsEnabled(self):
        r"""StaticCondensationIsEnabled(BilinearForm self) -> bool"""
        return _bilinearform.BilinearForm_StaticCondensationIsEnabled(self)

    def SCFESpace(self):
        r"""SCFESpace(BilinearForm self) -> FiniteElementSpace"""
        return _bilinearform.BilinearForm_SCFESpace(self)

    def EnableHybridization(self, constr_space, constr_integ, ess_tdof_list):
        r"""EnableHybridization(BilinearForm self, FiniteElementSpace constr_space, BilinearFormIntegrator constr_integ, intArray ess_tdof_list)"""
        val = _bilinearform.BilinearForm_EnableHybridization(self, constr_space, constr_integ, ess_tdof_list)

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(constr_integ)
        constr_integ.thisown = 0


        return val


    def UsePrecomputedSparsity(self, ps=1):
        r"""UsePrecomputedSparsity(BilinearForm self, int ps=1)"""
        return _bilinearform.BilinearForm_UsePrecomputedSparsity(self, ps)

    def UseSparsity(self, *args):
        r"""
        UseSparsity(BilinearForm self, int * I, int * J, bool isSorted)
        UseSparsity(BilinearForm self, SparseMatrix A)
        """
        return _bilinearform.BilinearForm_UseSparsity(self, *args)

    def AllocateMatrix(self):
        r"""AllocateMatrix(BilinearForm self)"""
        return _bilinearform.BilinearForm_AllocateMatrix(self)

    def GetDBFI(self):
        r"""GetDBFI(BilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.BilinearForm_GetDBFI(self)

    def GetBBFI(self):
        r"""GetBBFI(BilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.BilinearForm_GetBBFI(self)

    def GetBBFI_Marker(self):
        r"""GetBBFI_Marker(BilinearForm self) -> mfem::Array< mfem::Array< int > * > *"""
        return _bilinearform.BilinearForm_GetBBFI_Marker(self)

    def GetFBFI(self):
        r"""GetFBFI(BilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.BilinearForm_GetFBFI(self)

    def GetBFBFI(self):
        r"""GetBFBFI(BilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.BilinearForm_GetBFBFI(self)

    def GetBFBFI_Marker(self):
        r"""GetBFBFI_Marker(BilinearForm self) -> mfem::Array< mfem::Array< int > * > *"""
        return _bilinearform.BilinearForm_GetBFBFI_Marker(self)

    def __call__(self, i, j):
        r"""__call__(BilinearForm self, int i, int j) -> double const &"""
        return _bilinearform.BilinearForm___call__(self, i, j)

    def Elem(self, *args):
        r"""
        Elem(BilinearForm self, int i, int j) -> double
        Elem(BilinearForm self, int i, int j) -> double const &
        """
        return _bilinearform.BilinearForm_Elem(self, *args)

    def Mult(self, x, y):
        r"""Mult(BilinearForm self, Vector x, Vector y)"""
        return _bilinearform.BilinearForm_Mult(self, x, y)

    def FullMult(self, x, y):
        r"""FullMult(BilinearForm self, Vector x, Vector y)"""
        return _bilinearform.BilinearForm_FullMult(self, x, y)

    def AddMult(self, x, y, a=1.0):
        r"""AddMult(BilinearForm self, Vector x, Vector y, double const a=1.0)"""
        return _bilinearform.BilinearForm_AddMult(self, x, y, a)

    def FullAddMult(self, x, y):
        r"""FullAddMult(BilinearForm self, Vector x, Vector y)"""
        return _bilinearform.BilinearForm_FullAddMult(self, x, y)

    def AddMultTranspose(self, x, y, a=1.0):
        r"""AddMultTranspose(BilinearForm self, Vector x, Vector y, double const a=1.0)"""
        return _bilinearform.BilinearForm_AddMultTranspose(self, x, y, a)

    def FullAddMultTranspose(self, x, y):
        r"""FullAddMultTranspose(BilinearForm self, Vector x, Vector y)"""
        return _bilinearform.BilinearForm_FullAddMultTranspose(self, x, y)

    def MultTranspose(self, x, y):
        r"""MultTranspose(BilinearForm self, Vector x, Vector y)"""
        return _bilinearform.BilinearForm_MultTranspose(self, x, y)

    def InnerProduct(self, x, y):
        r"""InnerProduct(BilinearForm self, Vector x, Vector y) -> double"""
        return _bilinearform.BilinearForm_InnerProduct(self, x, y)

    def Inverse(self):
        r"""Inverse(BilinearForm self) -> MatrixInverse"""
        return _bilinearform.BilinearForm_Inverse(self)

    def Finalize(self, skip_zeros=1):
        r"""Finalize(BilinearForm self, int skip_zeros=1)"""
        return _bilinearform.BilinearForm_Finalize(self, skip_zeros)

    def SpMat(self, *args):
        r"""
        SpMat(BilinearForm self) -> SparseMatrix
        SpMat(BilinearForm self) -> SparseMatrix
        """
        val = _bilinearform.BilinearForm_SpMat(self, *args)

        if not hasattr(self, "_spmat"): self._spmat = []
        self._spmat.append(val)
        val.thisown=0 


        return val


    def LoseMat(self):
        r"""LoseMat(BilinearForm self) -> SparseMatrix"""
        return _bilinearform.BilinearForm_LoseMat(self)

    def SpMatElim(self, *args):
        r"""
        SpMatElim(BilinearForm self) -> SparseMatrix
        SpMatElim(BilinearForm self) -> SparseMatrix
        """
        return _bilinearform.BilinearForm_SpMatElim(self, *args)

    def AddDomainIntegrator(self, bfi):
        r"""AddDomainIntegrator(BilinearForm self, BilinearFormIntegrator bfi)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.BilinearForm_AddDomainIntegrator(self, bfi)


    def AddBoundaryIntegrator(self, *args):
        r"""
        AddBoundaryIntegrator(BilinearForm self, BilinearFormIntegrator bfi)
        AddBoundaryIntegrator(BilinearForm self, BilinearFormIntegrator bfi, intArray bdr_marker)
        """

        if not hasattr(self, "_integrators"): self._integrators = []
        bfi = args[0]	     
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.BilinearForm_AddBoundaryIntegrator(self, *args)


    def AddInteriorFaceIntegrator(self, bfi):
        r"""AddInteriorFaceIntegrator(BilinearForm self, BilinearFormIntegrator bfi)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.BilinearForm_AddInteriorFaceIntegrator(self, bfi)


    def AddBdrFaceIntegrator(self, *args):
        r"""
        AddBdrFaceIntegrator(BilinearForm self, BilinearFormIntegrator bfi)
        AddBdrFaceIntegrator(BilinearForm self, BilinearFormIntegrator bfi, intArray bdr_marker)
        """

        if not hasattr(self, "_integrators"): self._integrators = []
        bfi = args[0]
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.BilinearForm_AddBdrFaceIntegrator(self, *args)


    def Assemble(self, skip_zeros=1):
        r"""Assemble(BilinearForm self, int skip_zeros=1)"""
        return _bilinearform.BilinearForm_Assemble(self, skip_zeros)

    def GetProlongation(self):
        r"""GetProlongation(BilinearForm self) -> Operator"""
        return _bilinearform.BilinearForm_GetProlongation(self)

    def GetRestriction(self):
        r"""GetRestriction(BilinearForm self) -> Operator"""
        return _bilinearform.BilinearForm_GetRestriction(self)

    def RecoverFEMSolution(self, X, b, x):
        r"""RecoverFEMSolution(BilinearForm self, Vector X, Vector b, Vector x)"""
        return _bilinearform.BilinearForm_RecoverFEMSolution(self, X, b, x)

    def ComputeElementMatrices(self):
        r"""ComputeElementMatrices(BilinearForm self)"""
        return _bilinearform.BilinearForm_ComputeElementMatrices(self)

    def FreeElementMatrices(self):
        r"""FreeElementMatrices(BilinearForm self)"""
        return _bilinearform.BilinearForm_FreeElementMatrices(self)

    def ComputeElementMatrix(self, i, elmat):
        r"""ComputeElementMatrix(BilinearForm self, int i, DenseMatrix elmat)"""
        return _bilinearform.BilinearForm_ComputeElementMatrix(self, i, elmat)

    def AssembleElementMatrix(self, i, elmat, vdofs, skip_zeros=1):
        r"""AssembleElementMatrix(BilinearForm self, int i, DenseMatrix elmat, intArray vdofs, int skip_zeros=1)"""
        return _bilinearform.BilinearForm_AssembleElementMatrix(self, i, elmat, vdofs, skip_zeros)

    def AssembleBdrElementMatrix(self, i, elmat, vdofs, skip_zeros=1):
        r"""AssembleBdrElementMatrix(BilinearForm self, int i, DenseMatrix elmat, intArray vdofs, int skip_zeros=1)"""
        return _bilinearform.BilinearForm_AssembleBdrElementMatrix(self, i, elmat, vdofs, skip_zeros)

    def EliminateEssentialBC(self, *args):
        r"""
        EliminateEssentialBC(BilinearForm self, intArray bdr_attr_is_ess, Vector sol, Vector rhs, mfem::Matrix::DiagonalPolicy dpolicy=DIAG_ONE)
        EliminateEssentialBC(BilinearForm self, intArray bdr_attr_is_ess, mfem::Matrix::DiagonalPolicy dpolicy=DIAG_ONE)
        """
        return _bilinearform.BilinearForm_EliminateEssentialBC(self, *args)

    def EliminateEssentialBCDiag(self, bdr_attr_is_ess, value):
        r"""EliminateEssentialBCDiag(BilinearForm self, intArray bdr_attr_is_ess, double value)"""
        return _bilinearform.BilinearForm_EliminateEssentialBCDiag(self, bdr_attr_is_ess, value)

    def EliminateVDofs(self, *args):
        r"""
        EliminateVDofs(BilinearForm self, intArray vdofs, Vector sol, Vector rhs, mfem::Matrix::DiagonalPolicy dpolicy=DIAG_ONE)
        EliminateVDofs(BilinearForm self, intArray vdofs, mfem::Matrix::DiagonalPolicy dpolicy=DIAG_ONE)
        """
        return _bilinearform.BilinearForm_EliminateVDofs(self, *args)

    def EliminateEssentialBCFromDofs(self, *args):
        r"""
        EliminateEssentialBCFromDofs(BilinearForm self, intArray ess_dofs, Vector sol, Vector rhs, mfem::Matrix::DiagonalPolicy dpolicy=DIAG_ONE)
        EliminateEssentialBCFromDofs(BilinearForm self, intArray ess_dofs, mfem::Matrix::DiagonalPolicy dpolicy=DIAG_ONE)
        """
        return _bilinearform.BilinearForm_EliminateEssentialBCFromDofs(self, *args)

    def EliminateEssentialBCFromDofsDiag(self, ess_dofs, value):
        r"""EliminateEssentialBCFromDofsDiag(BilinearForm self, intArray ess_dofs, double value)"""
        return _bilinearform.BilinearForm_EliminateEssentialBCFromDofsDiag(self, ess_dofs, value)

    def EliminateVDofsInRHS(self, vdofs, x, b):
        r"""EliminateVDofsInRHS(BilinearForm self, intArray vdofs, Vector x, Vector b)"""
        return _bilinearform.BilinearForm_EliminateVDofsInRHS(self, vdofs, x, b)

    def FullInnerProduct(self, x, y):
        r"""FullInnerProduct(BilinearForm self, Vector x, Vector y) -> double"""
        return _bilinearform.BilinearForm_FullInnerProduct(self, x, y)

    def Update(self, nfes=None):
        r"""Update(BilinearForm self, FiniteElementSpace nfes=None)"""
        return _bilinearform.BilinearForm_Update(self, nfes)

    def GetFES(self):
        r"""GetFES(BilinearForm self) -> FiniteElementSpace"""
        return _bilinearform.BilinearForm_GetFES(self)

    def FESpace(self, *args):
        r"""
        FESpace(BilinearForm self) -> FiniteElementSpace
        FESpace(BilinearForm self) -> FiniteElementSpace
        """
        return _bilinearform.BilinearForm_FESpace(self, *args)

    def SetDiagonalPolicy(self, policy):
        r"""SetDiagonalPolicy(BilinearForm self, mfem::Matrix::DiagonalPolicy policy)"""
        return _bilinearform.BilinearForm_SetDiagonalPolicy(self, policy)
    __swig_destroy__ = _bilinearform.delete_BilinearForm

    def FormLinearSystem(self, *args):
        r"""
        FormLinearSystem(BilinearForm self, intArray ess_tdof_list, Vector x, Vector b, OperatorHandle A, Vector X, Vector B, int copy_interior=0)
        FormLinearSystem(BilinearForm self, intArray ess_tdof_list, Vector x, Vector b, SparseMatrix A, Vector X, Vector B, int copy_interior=0)
        """
        return _bilinearform.BilinearForm_FormLinearSystem(self, *args)

    def FormSystemMatrix(self, *args):
        r"""
        FormSystemMatrix(BilinearForm self, intArray ess_tdof_list, OperatorHandle A)
        FormSystemMatrix(BilinearForm self, intArray ess_tdof_list, SparseMatrix A)
        """
        return _bilinearform.BilinearForm_FormSystemMatrix(self, *args)
    def __disown__(self):
        self.this.disown()
        _bilinearform.disown_BilinearForm(self)
        return weakref.proxy(self)

# Register BilinearForm in _bilinearform:
_bilinearform.BilinearForm_swigregister(BilinearForm)

class MixedBilinearForm(mfem._ser.matrix.Matrix):
    r"""Proxy of C++ mfem::MixedBilinearForm class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(MixedBilinearForm self, FiniteElementSpace tr_fes, FiniteElementSpace te_fes) -> MixedBilinearForm
        __init__(MixedBilinearForm self, FiniteElementSpace tr_fes, FiniteElementSpace te_fes, MixedBilinearForm mbf) -> MixedBilinearForm
        """
        _bilinearform.MixedBilinearForm_swiginit(self, _bilinearform.new_MixedBilinearForm(*args))

    def Elem(self, *args):
        r"""
        Elem(MixedBilinearForm self, int i, int j) -> double
        Elem(MixedBilinearForm self, int i, int j) -> double const &
        """
        return _bilinearform.MixedBilinearForm_Elem(self, *args)

    def Mult(self, x, y):
        r"""Mult(MixedBilinearForm self, Vector x, Vector y)"""
        return _bilinearform.MixedBilinearForm_Mult(self, x, y)

    def AddMult(self, x, y, a=1.0):
        r"""AddMult(MixedBilinearForm self, Vector x, Vector y, double const a=1.0)"""
        return _bilinearform.MixedBilinearForm_AddMult(self, x, y, a)

    def AddMultTranspose(self, x, y, a=1.0):
        r"""AddMultTranspose(MixedBilinearForm self, Vector x, Vector y, double const a=1.0)"""
        return _bilinearform.MixedBilinearForm_AddMultTranspose(self, x, y, a)

    def MultTranspose(self, x, y):
        r"""MultTranspose(MixedBilinearForm self, Vector x, Vector y)"""
        return _bilinearform.MixedBilinearForm_MultTranspose(self, x, y)

    def Inverse(self):
        r"""Inverse(MixedBilinearForm self) -> MatrixInverse"""
        return _bilinearform.MixedBilinearForm_Inverse(self)

    def Finalize(self, skip_zeros=1):
        r"""Finalize(MixedBilinearForm self, int skip_zeros=1)"""
        return _bilinearform.MixedBilinearForm_Finalize(self, skip_zeros)

    def GetBlocks(self, blocks):
        r"""GetBlocks(MixedBilinearForm self, mfem::Array2D< mfem::SparseMatrix * > & blocks)"""
        return _bilinearform.MixedBilinearForm_GetBlocks(self, blocks)

    def SpMat(self, *args):
        r"""
        SpMat(MixedBilinearForm self) -> SparseMatrix
        SpMat(MixedBilinearForm self) -> SparseMatrix
        """
        val = _bilinearform.MixedBilinearForm_SpMat(self, *args)

        if not hasattr(self, "_spmat"): self._spmat = []
        self._spmat.append(val)
        val.thisown=0 


        return val


    def LoseMat(self):
        r"""LoseMat(MixedBilinearForm self) -> SparseMatrix"""
        return _bilinearform.MixedBilinearForm_LoseMat(self)

    def AddDomainIntegrator(self, bfi):
        r"""AddDomainIntegrator(MixedBilinearForm self, BilinearFormIntegrator bfi)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.MixedBilinearForm_AddDomainIntegrator(self, bfi)


    def AddBoundaryIntegrator(self, bfi):
        r"""AddBoundaryIntegrator(MixedBilinearForm self, BilinearFormIntegrator bfi)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.MixedBilinearForm_AddBoundaryIntegrator(self, bfi)


    def AddTraceFaceIntegrator(self, bfi):
        r"""AddTraceFaceIntegrator(MixedBilinearForm self, BilinearFormIntegrator bfi)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(bfi)
        bfi.thisown=0 


        return _bilinearform.MixedBilinearForm_AddTraceFaceIntegrator(self, bfi)


    def GetDBFI(self):
        r"""GetDBFI(MixedBilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.MixedBilinearForm_GetDBFI(self)

    def GetBBFI(self):
        r"""GetBBFI(MixedBilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.MixedBilinearForm_GetBBFI(self)

    def GetTFBFI(self):
        r"""GetTFBFI(MixedBilinearForm self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.MixedBilinearForm_GetTFBFI(self)

    def Assemble(self, skip_zeros=1):
        r"""Assemble(MixedBilinearForm self, int skip_zeros=1)"""
        return _bilinearform.MixedBilinearForm_Assemble(self, skip_zeros)

    def ConformingAssemble(self):
        r"""ConformingAssemble(MixedBilinearForm self)"""
        return _bilinearform.MixedBilinearForm_ConformingAssemble(self)

    def EliminateTrialDofs(self, bdr_attr_is_ess, sol, rhs):
        r"""EliminateTrialDofs(MixedBilinearForm self, intArray bdr_attr_is_ess, Vector sol, Vector rhs)"""
        return _bilinearform.MixedBilinearForm_EliminateTrialDofs(self, bdr_attr_is_ess, sol, rhs)

    def EliminateEssentialBCFromTrialDofs(self, marked_vdofs, sol, rhs):
        r"""EliminateEssentialBCFromTrialDofs(MixedBilinearForm self, intArray marked_vdofs, Vector sol, Vector rhs)"""
        return _bilinearform.MixedBilinearForm_EliminateEssentialBCFromTrialDofs(self, marked_vdofs, sol, rhs)

    def EliminateTestDofs(self, bdr_attr_is_ess):
        r"""EliminateTestDofs(MixedBilinearForm self, intArray bdr_attr_is_ess)"""
        return _bilinearform.MixedBilinearForm_EliminateTestDofs(self, bdr_attr_is_ess)

    def Update(self):
        r"""Update(MixedBilinearForm self)"""
        return _bilinearform.MixedBilinearForm_Update(self)
    __swig_destroy__ = _bilinearform.delete_MixedBilinearForm

# Register MixedBilinearForm in _bilinearform:
_bilinearform.MixedBilinearForm_swigregister(MixedBilinearForm)

class DiscreteLinearOperator(MixedBilinearForm):
    r"""Proxy of C++ mfem::DiscreteLinearOperator class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, domain_fes, range_fes):
        r"""__init__(DiscreteLinearOperator self, FiniteElementSpace domain_fes, FiniteElementSpace range_fes) -> DiscreteLinearOperator"""
        _bilinearform.DiscreteLinearOperator_swiginit(self, _bilinearform.new_DiscreteLinearOperator(domain_fes, range_fes))

    def AddDomainInterpolator(self, di):
        r"""AddDomainInterpolator(DiscreteLinearOperator self, DiscreteInterpolator di)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(di)
        di.thisown=0 


        return _bilinearform.DiscreteLinearOperator_AddDomainInterpolator(self, di)


    def AddTraceFaceInterpolator(self, di):
        r"""AddTraceFaceInterpolator(DiscreteLinearOperator self, DiscreteInterpolator di)"""

        if not hasattr(self, "_integrators"): self._integrators = []
        self._integrators.append(di)
        di.thisown=0 


        return _bilinearform.DiscreteLinearOperator_AddTraceFaceInterpolator(self, di)


    def GetDI(self):
        r"""GetDI(DiscreteLinearOperator self) -> mfem::Array< mfem::BilinearFormIntegrator * > *"""
        return _bilinearform.DiscreteLinearOperator_GetDI(self)

    def Assemble(self, skip_zeros=1):
        r"""Assemble(DiscreteLinearOperator self, int skip_zeros=1)"""
        return _bilinearform.DiscreteLinearOperator_Assemble(self, skip_zeros)
    __swig_destroy__ = _bilinearform.delete_DiscreteLinearOperator

# Register DiscreteLinearOperator in _bilinearform:
_bilinearform.DiscreteLinearOperator_swigregister(DiscreteLinearOperator)



