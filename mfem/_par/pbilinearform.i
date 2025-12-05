//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") pbilinearform
%{
//#include <mpi.h>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "handle.i"
%import "bilinearform.i"
%import "pfespace.i"
%import "pgridfunc.i"
%import "hypre.i"
%import "../common/exception.i"

%newobject mfem::ParBilinearForm::ParallelAssemble;

%include "fem/pbilinearform.hpp"

// This template is defined in baseclass.
// We need to add it to make template macro works.
%extend mfem::ParBilinearForm {
   template <typename OpType>
     void FormLinearSystem(const mfem::Array<int> &ess_tdof_list,
			   mfem::Vector &x, mfem::Vector &b,
                           OpType &A, mfem::Vector &X, mfem::Vector &B,
                           int copy_interior = 0){
     return self->mfem::BilinearForm::FormLinearSystem(ess_tdof_list,
						       x, b, A, X, B,
						       copy_interior);
   }
   template <typename OpType>
     void FormSystemMatrix(const mfem::Array<int> &ess_tdof_list,
			   OpType &A){
     return self->mfem::BilinearForm::FormSystemMatrix(ess_tdof_list, A);
   }
};
%extend mfem::ParMixedBilinearForm {
   template <typename OpType>
     void FormRectangularSystemMatrix(const Array<int> &trial_tdof_list,
				      const Array<int> &test_tdof_list, OpType &A){
     return self->mfem::MixedBilinearForm::FormRectangularSystemMatrix(trial_tdof_list,
								       test_tdof_list,
								       A);
   }

 template <typename OpType>
   void FormRectangularLinearSystem(const Array<int> &trial_tdof_list,
                                    const Array<int> &test_tdof_list,
                                    Vector &x, Vector &b,
                                    OpType &A, Vector &X, Vector &B){
     return self->mfem::MixedBilinearForm::FormRectangularLinearSystem(trial_tdof_list,
								       test_tdof_list,
								       x, b,
								       A,
								       X, B);
 }
};

// instatitate template methods

%ignore FORM_SYSTEM_MATRIX_WRAP;

%define P_FORM_SYSTEM_MATRIX_WRAP(OsType)
%template(FormLinearSystem) mfem::ParBilinearForm::FormLinearSystem<OsType>;
%template(FormSystemMatrix) mfem::ParBilinearForm::FormSystemMatrix<OsType>;
%template(FormRectangularLinearSystem) mfem::ParMixedBilinearForm::FormRectangularLinearSystem<OsType>;
%template(FormRectangularSystemMatrix) mfem::ParMixedBilinearForm::FormRectangularSystemMatrix<OsType>;
%enddef

P_FORM_SYSTEM_MATRIX_WRAP(mfem::SparseMatrix)

#ifdef MFEM_USE_MPI
  P_FORM_SYSTEM_MATRIX_WRAP(mfem::HypreParMatrix)
#endif

#ifdef MFEM_USE_PETSC
  P_FORM_SYSTEM_MATRIX_WRAP(mfem::PetscParMatrix)
#endif

