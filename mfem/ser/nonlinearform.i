%module nonlinearform
%{
#include "fem/linearform.hpp"    
#include "fem/nonlininteg.hpp"
#include "fem/nonlinearform.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"     
%}

%init %{
import_array();
%}

%include "exception.i"
%import "operators.i"
%import "fespace.i"
%import "nonlininteg.i"
%import "../common/exception.i"
%include "../common/typemap_macros.i"


LIST_TO_OBJARRAY_IN(mfem::Array<mfem::FiniteElementSpace *> &f,
		    mfem::FiniteElementSpace *)
LIST_TO_OBJARRAY_IN(const mfem::Array<mfem::Array<int> *> &bdr_attr_is_ess,
		    mfem::Array<int> *)
LIST_TO_OBJARRAY_IN(mfem::Array<mfem::Vector *> &rhs,
		    mfem::Vector *)

namespace mfem { 
%pythonprepend NonlinearForm::AddDomainIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
#    self._integrators.append(nlfi)
    nlfi.thisown=0 
%}
%pythonprepend NonlinearForm::AddInteriorFaceIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
#    self._integrators.append(nlfi)
    nlfi.thisown=0 
%}
%pythonprepend NonlinearForm::AddBdrFaceIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
#    self._integrators.append(nlfi)
    nlfi = args[0]
    nlfi.thisown=0 
%}
}


%include "fem/nonlinearform.hpp"

%extend mfem::NonlinearForm{
  mfem::SparseMatrix *GetGradientMatrix(const mfem::Vector &x) const
  {
      return  dynamic_cast<mfem::SparseMatrix *>(&$self->GetGradient(x));
  }
};
