%module(package="mfem._par") pfespace
%{
#include <mpi.h>
#include  "config/config.hpp"
#include "fem/linearform.hpp"  
#include "fem/pfespace.hpp"
#include "fem/restriction.hpp"
#include "fem/prestriction.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"    
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array();
%}

%include "exception.i"
 //%include "../common/cpointers.i"
 //%import "cpointers.i"
%import "operators.i"
%import "fespace.i"
%import "pmesh.i"
%import "hypre.i"
%import "restriction.i"
%import "../common/exception.i"

%import "../common/hypre_int.i"

%pointer_class(int, intp);
%immutable face_nbr_glob_dof_map;
//DoF accesser

// default number is -1, which conflict with error code of PyArray_PyIntAsInt...
%typemap(typecheck) (int component = -1) {
   $1 = PyInt_Check($input) ? 1 : 0;
}

%feature("shadow") mfem::ParFiniteElementSpace::GetSharedEdgeDofs %{
def GetSharedEdgeDofs(self, group, ei):
    from  .array import intArray
    dofs = intArray() 
    $action(self, group, ei, dofs)
    return dofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetSharedFaceDofs %{
def GetSharedFaceDofs(self, group, fi):
    from  .array import intArray
    dofs = intArray()
    $action(self, group, fi, dofs)
    return dofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetSharedTriangleDofs %{
def GetSharedTriangleDofs(self, group, fi):
    from  .array import intArray
    dofs = intArray()
    $action(self, group, fi, dofs)
    return dofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetSharedQuadrilateralDofs %{
def GetSharedQuadrilateralDofs(self, group, fi):
    from  .array import intArray
    dofs = intArray()
    $action(self, group, fi, dofs)
    return dofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetFaceNbrElementVDofs %{
def GetFaceNbrElementVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    $action(self, i, vdofs)
    return vdofs.ToList()
%}

%include "fem/pfespace.hpp"

%extend mfem::ParFiniteElementSpace{
  virtual DofTransformation *GetElementDofTransformation(int elem) const{
    mfem::Array<int> dofs;
    return self->GetElementDofs(elem, dofs);
  }
  virtual DofTransformation *GetBdrElementDofTransformation(int bel) const {
    mfem::Array<int> dofs;
    return self->GetBdrElementDofs(bel, dofs);
  }
  virtual DofTransformation *GetFaceNbrVDofTransformation(int elem) const {
    mfem::Array<int> dofs;    
    return self->GetFaceNbrElementVDofs(elem, dofs);
  }
};

