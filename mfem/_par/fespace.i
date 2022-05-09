%module(package="mfem._par") fespace

%feature("autodoc", "1");

%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "fem/fem.hpp"
#include "fem/fe_coll.hpp"
#include "fem/fespace.hpp"
#include "fem/eltrans.hpp"
#include "fem/coefficient.hpp"
#include "fem/intrules.hpp"  
#include "fem/restriction.hpp"
#include "../common/io_stream.hpp"      
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"        
%}

%init %{
import_array();
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "coefficient.i"
%import "matrix.i"
%import "mesh.i"
%import "intrules.i"
%import "fe.i"
%import "fe_coll.i"
%import "doftrans.i"
%import "densemat.i"
%import "sparsemat.i"
%import "eltrans.i"
%import "lininteg.i"
%import "handle.i"
%import "restriction.i"
%include "../common/typemap_macros.i"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

// default number is -1, which conflict with error code of PyArray_PyIntAsInt...
INT_DEFAULT_NEGATIVE_ONE(int ndofs = -1)
INT_DEFAULT_NEGATIVE_ONE(int component = -1)

//VDoF accesser
%feature("shadow") mfem::FiniteElementSpace::GetBdrElementVDofs %{
def GetBdrElementVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementVDofs %{
def GetElementVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetFaceVDofs %{
def GetFaceVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeVDofs %{
def GetEdgeVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetVertexVDofs %{
def GetVertexVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementInteriorVDofs %{
def GetElementInteriorVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeInteriorVDofs %{
def GetEdgeInteriorVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}

//DoF accesser
%feature("shadow") mfem::FiniteElementSpace::GetBdrElementDofs %{
def GetBdrElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementDofs %{
def GetElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetFaceDofs %{
def GetFaceDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeDofs %{
def GetEdgeDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetVertexDofs %{
def GetVertexDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementInteriorDofs %{
def GetElementInteriorDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeInteriorDofs %{
def GetEdgeInteriorDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}

%pythonappend mfem::FiniteElementSpace::FiniteElementSpace%{

    '''
        FiniteElementSpace(Mesh *m, const FiniteElementCollection *f,
                      int vdim = 1, int ordering = Ordering::byNODES);
        keep reference to mesh and fec to prevent it from
        freeed from pytho garbage collector
    '''
    if hasattr(self, 'this'):
        self.mesh = args[0]
        self.fec = args[1]
      
%}

/* define FiniteElementSpacePtrArray */
%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::FiniteElementSpace *, 1)

%import "../common/data_size_typemap.i"
XXXPTR_SIZE_IN(mfem::FiniteElementSpace **data_, int asize, mfem::FiniteElementSpace *)

%import "../common/array_instantiation_macro.i"
IGNORE_ARRAY_METHODS(mfem::FiniteElementSpace *)
INSTANTIATE_ARRAY0(FiniteElementSpace *, FiniteElementSpace, 1)

%include "fem/fespace.hpp"

/*
fem/fespace.hpp:   void Save(std::ostream &out) const;
fem/fespace.hpp:   void Save(std::ostream &out) const;
*/
OSTREAM_ADD_DEFAULT_STDOUT_FILE(FiniteElementSpace, Save)
OSTREAM_ADD_DEFAULT_STDOUT_FILE(QuadratureSpace, Save)

%extend mfem::FiniteElementSpace{
  virtual DofTransformation *GetElementDofTransformation(int elem) const{
    mfem::Array<int> dofs;
    return self->GetElementDofs(elem, dofs);
  }
  virtual DofTransformation *GetBdrElementDofTransformation(int bel) const {
    mfem::Array<int> dofs;
    return self->GetBdrElementDofs(bel, dofs);
  }
  virtual DofTransformation *GetElementVDofTransformation(int elem) const {
    mfem::Array<int> dofs;    
    return self->GetElementVDofs(elem, dofs);
  }
  virtual DofTransformation *GetBdrElementVDofTransformation(int bel) const {
    mfem::Array<int> dofs;        
    return self->GetBdrElementVDofs(bel, dofs);
  }
};
