%module  handle
%{
#include "config/config.hpp"  
#include "linalg/hypre.hpp"
#include "linalg/handle.hpp"  
#include "fem/gridfunc.hpp"  
#include "fem/linearform.hpp"
#include "fem/pfespace.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
%}
%include  "config/_config.hpp" // include mfem MACRO

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}

%import operators.i
#ifdef MFEM_USE_MPI
%import hypre.i
#endif
#ifdef MFEM_USE_PETSC
%include petsc.i
#endif

//
%pointer_class(int, intp);


%include "linalg/handle.hpp"

// instatitate template methods (step 1: Macro definition)
%define OPERATORHANDLE_WRAP(T)
%template(OperatorHandle) mfem::OperatorHandle::OperatorHandle<T >;
%enddef
%define AS_WRAP(T)
%template(As) mfem::OperatorHandle::As<T >;
%enddef
%define IS_WRAP(T)
%template(Is) mfem::OperatorHandle::Is<T >;
%enddef
%define GET_WRAP(T)
%template(Get) mfem::OperatorHandle::Get<T >;
%enddef
%define RESET_WRAP(T)
%template(Reset) mfem::OperatorHandle::Reset<T >;
%enddef
%define CONVERT_FROM_WRAP(T)
%template(ConvertFrom) mfem::OperatorHandle::ConvertFrom<T >;
%enddef

// instatitate template methods (step 2: Instantiation)
AS_WRAP(mfem::HypreParMatrix)
IS_WRAP(mfem::HypreParMatrix)
GET_WRAP(mfem::HypreParMatrix)
RESET_WRAP(mfem::HypreParMatrix)          
CONVERT_FROM_WRAP(mfem::HypreParMatrix)

#ifdef MFEM_USE_PETSC
AS_WRAP(mfem::PetscParMatrix)
IS_WRAP(mfem::PetscParMatrix)
GET_WRAP(mfem::PetscParMatrix)
RESET_WRAP(mfem::PetscParMatrix)          
CONVERT_FROM_WRAP(mfem::PetscParMatrix)
#endif     
