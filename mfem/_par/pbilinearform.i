%module(package="mfem._par") pbilinearform
%{
  //#include <mpi.h>  
#include "config/config.hpp"    
#include "fem/pbilinearform.hpp"
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
%import "handle.i"
%import "bilinearform.i"
%import "pfespace.i"
%import "hypre.i"
%import "../common/exception.i"

%newobject mfem::ParBilinearForm::ParallelAssemble;
%pointer_class(int, intp);

%include "fem/pbilinearform.hpp"

// This template is defined in baseclass.
// We need to add it to make template macro works.
%extend mfem::ParBilinearForm {
   template <typename OpType>
     void FormLinearSystem(const mfem::Array<int> &ess_tdof_list,
			   mfem::Vector &x, mfem::Vector &b,
                           OpType &A, mfem::Vector &X, mfem::Vector &B,
                           int copy_interior = 0){
     return self->mfem::BilinearForm::FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior);
   }
   template <typename OpType>
     void FormSystemMatrix(const mfem::Array<int> &ess_tdof_list,
			   OpType &A){
     return self->mfem::BilinearForm::FormSystemMatrix(ess_tdof_list, A);
   }
};  

// instatitate template methods

%ignore FORM_SYSTEM_MATRIX_WRAP;

%define P_FORM_SYSTEM_MATRIX_WRAP(OsType)
%template(FormLinearSystem) mfem::ParBilinearForm::FormLinearSystem<OsType>;
%template(FormSystemMatrix) mfem::ParBilinearForm::FormSystemMatrix<OsType>;
%enddef

P_FORM_SYSTEM_MATRIX_WRAP(mfem::SparseMatrix)
  
#ifdef MFEM_USE_MPI
  P_FORM_SYSTEM_MATRIX_WRAP(mfem::HypreParMatrix)
#endif
  
#ifdef MFEM_USE_PETSC
  P_FORM_SYSTEM_MATRIX_WRAP(mfem::PetscParMatrix)
#endif  

