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
    from . import _ncmesh
else:
    import _ncmesh

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
class Refinement(object):
    r"""Proxy of C++ mfem::Refinement class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    index = property(_ncmesh.Refinement_index_get, _ncmesh.Refinement_index_set, doc=r"""index : int""")
    ref_type = property(_ncmesh.Refinement_ref_type_get, _ncmesh.Refinement_ref_type_set, doc=r"""ref_type : char""")

    def __init__(self, *args):
        r"""
        __init__(Refinement self) -> Refinement
        __init__(Refinement self, int index, int type=7) -> Refinement
        """
        _ncmesh.Refinement_swiginit(self, _ncmesh.new_Refinement(*args))
    __swig_destroy__ = _ncmesh.delete_Refinement

# Register Refinement in _ncmesh:
_ncmesh.Refinement_swigregister(Refinement)

class Embedding(object):
    r"""Proxy of C++ mfem::Embedding class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    parent = property(_ncmesh.Embedding_parent_get, _ncmesh.Embedding_parent_set, doc=r"""parent : int""")
    matrix = property(_ncmesh.Embedding_matrix_get, _ncmesh.Embedding_matrix_set, doc=r"""matrix : int""")

    def __init__(self, *args):
        r"""
        __init__(Embedding self) -> Embedding
        __init__(Embedding self, int elem, int matrix=0) -> Embedding
        """
        _ncmesh.Embedding_swiginit(self, _ncmesh.new_Embedding(*args))
    __swig_destroy__ = _ncmesh.delete_Embedding

# Register Embedding in _ncmesh:
_ncmesh.Embedding_swigregister(Embedding)

class CoarseFineTransformations(object):
    r"""Proxy of C++ mfem::CoarseFineTransformations class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    point_matrices = property(_ncmesh.CoarseFineTransformations_point_matrices_get, _ncmesh.CoarseFineTransformations_point_matrices_set, doc=r"""point_matrices : std::map<(mfem::Geometry::Type,mfem::DenseTensor)>""")
    embeddings = property(_ncmesh.CoarseFineTransformations_embeddings_get, doc=r"""embeddings : mfem::Array<(mfem::Embedding)>""")

    def GetPointMatrices(self, geom):
        r"""GetPointMatrices(CoarseFineTransformations self, mfem::Geometry::Type geom) -> DenseTensor"""
        return _ncmesh.CoarseFineTransformations_GetPointMatrices(self, geom)

    def GetCoarseToFineMap(self, fine_mesh, coarse_to_fine, coarse_to_ref_type, ref_type_to_matrix, ref_type_to_geom):
        r"""GetCoarseToFineMap(CoarseFineTransformations self, Mesh fine_mesh, Table coarse_to_fine, intArray coarse_to_ref_type, Table ref_type_to_matrix, mfem::Array< mfem::Geometry::Type > & ref_type_to_geom)"""
        return _ncmesh.CoarseFineTransformations_GetCoarseToFineMap(self, fine_mesh, coarse_to_fine, coarse_to_ref_type, ref_type_to_matrix, ref_type_to_geom)

    def Clear(self):
        r"""Clear(CoarseFineTransformations self)"""
        return _ncmesh.CoarseFineTransformations_Clear(self)

    def MemoryUsage(self):
        r"""MemoryUsage(CoarseFineTransformations self) -> long"""
        return _ncmesh.CoarseFineTransformations_MemoryUsage(self)

    def __init__(self):
        r"""__init__(CoarseFineTransformations self) -> CoarseFineTransformations"""
        _ncmesh.CoarseFineTransformations_swiginit(self, _ncmesh.new_CoarseFineTransformations())
    __swig_destroy__ = _ncmesh.delete_CoarseFineTransformations

# Register CoarseFineTransformations in _ncmesh:
_ncmesh.CoarseFineTransformations_swigregister(CoarseFineTransformations)

class NCMesh(object):
    r"""Proxy of C++ mfem::NCMesh class."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(NCMesh self, Mesh mesh, std::istream * vertex_parents=None) -> NCMesh
        __init__(NCMesh self, NCMesh other) -> NCMesh
        """
        _ncmesh.NCMesh_swiginit(self, _ncmesh.new_NCMesh(*args))
    __swig_destroy__ = _ncmesh.delete_NCMesh

    def Dimension(self):
        r"""Dimension(NCMesh self) -> int"""
        return _ncmesh.NCMesh_Dimension(self)

    def SpaceDimension(self):
        r"""SpaceDimension(NCMesh self) -> int"""
        return _ncmesh.NCMesh_SpaceDimension(self)

    def GetNVertices(self):
        r"""GetNVertices(NCMesh self) -> int"""
        return _ncmesh.NCMesh_GetNVertices(self)

    def GetNEdges(self):
        r"""GetNEdges(NCMesh self) -> int"""
        return _ncmesh.NCMesh_GetNEdges(self)

    def GetNFaces(self):
        r"""GetNFaces(NCMesh self) -> int"""
        return _ncmesh.NCMesh_GetNFaces(self)

    def Refine(self, refinements):
        r"""Refine(NCMesh self, mfem::Array< mfem::Refinement > const & refinements)"""
        return _ncmesh.NCMesh_Refine(self, refinements)

    def LimitNCLevel(self, max_nc_level):
        r"""LimitNCLevel(NCMesh self, int max_nc_level)"""
        return _ncmesh.NCMesh_LimitNCLevel(self, max_nc_level)

    def GetDerefinementTable(self):
        r"""GetDerefinementTable(NCMesh self) -> Table"""
        return _ncmesh.NCMesh_GetDerefinementTable(self)

    def CheckDerefinementNCLevel(self, deref_table, level_ok, max_nc_level):
        r"""CheckDerefinementNCLevel(NCMesh self, Table deref_table, intArray level_ok, int max_nc_level)"""
        return _ncmesh.NCMesh_CheckDerefinementNCLevel(self, deref_table, level_ok, max_nc_level)

    def Derefine(self, derefs):
        r"""Derefine(NCMesh self, intArray derefs)"""
        return _ncmesh.NCMesh_Derefine(self, derefs)

    def GetFaceList(self):
        r"""GetFaceList(NCMesh self) -> mfem::NCMesh::NCList const &"""
        return _ncmesh.NCMesh_GetFaceList(self)

    def GetEdgeList(self):
        r"""GetEdgeList(NCMesh self) -> mfem::NCMesh::NCList const &"""
        return _ncmesh.NCMesh_GetEdgeList(self)

    def GetVertexList(self):
        r"""GetVertexList(NCMesh self) -> mfem::NCMesh::NCList const &"""
        return _ncmesh.NCMesh_GetVertexList(self)

    def GetNCList(self, entity):
        r"""GetNCList(NCMesh self, int entity) -> mfem::NCMesh::NCList const &"""
        return _ncmesh.NCMesh_GetNCList(self, entity)

    def MarkCoarseLevel(self):
        r"""MarkCoarseLevel(NCMesh self)"""
        return _ncmesh.NCMesh_MarkCoarseLevel(self)

    def GetRefinementTransforms(self):
        r"""GetRefinementTransforms(NCMesh self) -> CoarseFineTransformations"""
        return _ncmesh.NCMesh_GetRefinementTransforms(self)

    def GetDerefinementTransforms(self):
        r"""GetDerefinementTransforms(NCMesh self) -> CoarseFineTransformations"""
        return _ncmesh.NCMesh_GetDerefinementTransforms(self)

    def ClearTransforms(self):
        r"""ClearTransforms(NCMesh self)"""
        return _ncmesh.NCMesh_ClearTransforms(self)

    @staticmethod
    def GridSfcOrdering2D(width, height, coords):
        r"""GridSfcOrdering2D(int width, int height, intArray coords)"""
        return _ncmesh.NCMesh_GridSfcOrdering2D(width, height, coords)

    @staticmethod
    def GridSfcOrdering3D(width, height, depth, coords):
        r"""GridSfcOrdering3D(int width, int height, int depth, intArray coords)"""
        return _ncmesh.NCMesh_GridSfcOrdering3D(width, height, depth, coords)

    def GetEdgeVertices(self, edge_id, vert_index, oriented=True):
        r"""GetEdgeVertices(NCMesh self, mfem::NCMesh::MeshId const & edge_id, int [2] vert_index, bool oriented=True)"""
        return _ncmesh.NCMesh_GetEdgeVertices(self, edge_id, vert_index, oriented)

    def GetEdgeNCOrientation(self, edge_id):
        r"""GetEdgeNCOrientation(NCMesh self, mfem::NCMesh::MeshId const & edge_id) -> int"""
        return _ncmesh.NCMesh_GetEdgeNCOrientation(self, edge_id)

    def GetFaceVerticesEdges(self, face_id, vert_index, edge_index, edge_orientation):
        r"""GetFaceVerticesEdges(NCMesh self, mfem::NCMesh::MeshId const & face_id, int [4] vert_index, int [4] edge_index, int [4] edge_orientation)"""
        return _ncmesh.NCMesh_GetFaceVerticesEdges(self, face_id, vert_index, edge_index, edge_orientation)

    def GetEdgeMaster(self, v1, v2):
        r"""GetEdgeMaster(NCMesh self, int v1, int v2) -> int"""
        return _ncmesh.NCMesh_GetEdgeMaster(self, v1, v2)

    def GetBoundaryClosure(self, bdr_attr_is_ess, bdr_vertices, bdr_edges):
        r"""GetBoundaryClosure(NCMesh self, intArray bdr_attr_is_ess, intArray bdr_vertices, intArray bdr_edges)"""
        return _ncmesh.NCMesh_GetBoundaryClosure(self, bdr_attr_is_ess, bdr_vertices, bdr_edges)

    def GetElementGeometry(self):
        r"""GetElementGeometry(NCMesh self) -> mfem::Geometry::Type"""
        return _ncmesh.NCMesh_GetElementGeometry(self)

    def GetFaceGeometry(self):
        r"""GetFaceGeometry(NCMesh self) -> mfem::Geometry::Type"""
        return _ncmesh.NCMesh_GetFaceGeometry(self)

    def GetElementDepth(self, i):
        r"""GetElementDepth(NCMesh self, int i) -> int"""
        return _ncmesh.NCMesh_GetElementDepth(self, i)

    def LoadVertexParents(self, input):
        r"""LoadVertexParents(NCMesh self, std::istream & input)"""
        return _ncmesh.NCMesh_LoadVertexParents(self, input)

    def LoadCoarseElements(self, input):
        r"""LoadCoarseElements(NCMesh self, std::istream & input)"""
        return _ncmesh.NCMesh_LoadCoarseElements(self, input)

    def SetVertexPositions(self, vertices):
        r"""SetVertexPositions(NCMesh self, mfem::Array< mfem::Vertex > const & vertices)"""
        return _ncmesh.NCMesh_SetVertexPositions(self, vertices)

    def Trim(self):
        r"""Trim(NCMesh self)"""
        return _ncmesh.NCMesh_Trim(self)

    def MemoryUsage(self):
        r"""MemoryUsage(NCMesh self) -> long"""
        return _ncmesh.NCMesh_MemoryUsage(self)

    def PrintMemoryDetail(self):
        r"""PrintMemoryDetail(NCMesh self) -> int"""
        return _ncmesh.NCMesh_PrintMemoryDetail(self)

    def PrintVertexParents(self, *args):
        r"""
        PrintVertexParents(NCMesh self, std::ostream & out)
        PrintVertexParents(NCMesh self, char const * file, int precision=8)
        PrintVertexParents(NCMesh self)
        """
        return _ncmesh.NCMesh_PrintVertexParents(self, *args)

    def PrintCoarseElements(self, *args):
        r"""
        PrintCoarseElements(NCMesh self, std::ostream & out)
        PrintCoarseElements(NCMesh self, char const * file, int precision=8)
        PrintCoarseElements(NCMesh self)
        """
        return _ncmesh.NCMesh_PrintCoarseElements(self, *args)

    def PrintStats(self, *args):
        r"""
        PrintStats(NCMesh self, std::ostream & out=mfem::out)
        PrintStats(NCMesh self, char const * file, int precision=8)
        """
        return _ncmesh.NCMesh_PrintStats(self, *args)

# Register NCMesh in _ncmesh:
_ncmesh.NCMesh_swigregister(NCMesh)

def NCMesh_GridSfcOrdering2D(width, height, coords):
    r"""NCMesh_GridSfcOrdering2D(int width, int height, intArray coords)"""
    return _ncmesh.NCMesh_GridSfcOrdering2D(width, height, coords)

def NCMesh_GridSfcOrdering3D(width, height, depth, coords):
    r"""NCMesh_GridSfcOrdering3D(int width, int height, int depth, intArray coords)"""
    return _ncmesh.NCMesh_GridSfcOrdering3D(width, height, depth, coords)



