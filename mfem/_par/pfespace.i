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

%feature("shadow") mfem::ParFiniteElementSpace::GetBdrElementDofs %{
def GetBdrElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _pfespace.ParFiniteElementSpace_GetBdrElementDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetElementDofs %{
def GetElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _pfespace.ParFiniteElementSpace_GetElementDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetFaceDofs %{
def GetFaceDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _pfespace.ParFiniteElementSpace_GetFaceDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetSharedEdgeDofs %{
  def GetSharedEdgeDofs(self, group, ei):
    from  .array import intArray
    dofs = intArray()      
    _pfespace.ParFiniteElementSpace_GetSharedEdgeDofs(self, group, ei, dofs)
    return dofs.ToList()      
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetSharedFaceDofs %{
  def GetSharedFaceDofs(self, group, fi):
    from  .array import intArray
    dofs = intArray()      
    _pfespace.ParFiniteElementSpace_GetSharedFaceDofs(self, group, fi, dofs)
    return dofs.ToList()      
%}

%include "fem/pfespace.hpp"
