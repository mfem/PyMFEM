%module(package="mfem._par") pbilinearform
%{
  //#include <mpi.h>  
#include "config/config.hpp"    
#include "fem/pbilinearform.hpp"
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
%import "handle.i"
%import "bilinearform.i"
%import "pfespace.i"
%import "hypre.i"
%import "../common/exception.i"

%newobject mfem::ParBilinearForm::ParallelAssemble;
%pointer_class(int, intp);

%include "fem/pbilinearform.hpp"

// instatitate template methods (step 1: Macro definition)
%define FORM_LINEAR_SYSTEM_WRAP(T)
%template(FormLinearSystem) mfem::ParBilinearForm::FormLinearSystem<T >;
%enddef
%define FORM_SYSTEM_MATRIX_WRAP(T)
%template(FormSystemMatrix) mfem::ParBilinearForm::FormSystemMatrix<T >;
%enddef

// instatitate template methods (step 2: Instantiation)
FORM_LINEAR_SYSTEM_WRAP(mfem::HypreParMatrix)
FORM_SYSTEM_MATRIX_WRAP(mfem::HypreParMatrix)
     
#ifdef MFEM_USE_PETSC
FORM_LINEAR_SYSTEM_WRAP(mfem::PetscParMatrix)
FORM_SYSTEM_MATRIX_WRAP(mfem::PetscParMatrix)
#endif          
 /*
(note above is the same as this)
%template(FormLinearSystem) mfem::ParBilinearForm::FormLinearSystem<mfem::HypreParMatrix>;
%template(FormSystemMatrix) mfem::ParBilinearForm::FormSystemMatrix<mfem::HypreParMatrix>;
 */
