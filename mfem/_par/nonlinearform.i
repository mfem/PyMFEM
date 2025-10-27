//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") nonlinearform
%{
#include <mpi.h>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
%}


%init %{
import_array1(-1);
%}

%include "exception.i"
%import "operators.i"
%import "fespace.i"
%import "bilinearform.i"
%import "nonlininteg.i"
%import "../common/exception.i"
%include "../common/typemap_macros.i"

namespace mfem {
%pythonprepend NonlinearForm::AddDomainIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
    nlfi = args[0]
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
    nlfi = args[0]
#    self._integrators.append(nlfi)
    nlfi.thisown=0
%}
%pythonprepend BlockNonlinearForm::AddDomainIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
    nlfi = args[0]
#    self._integrators.append(nlfi)
    nlfi.thisown=0
%}
%pythonprepend BlockNonlinearForm::AddInteriorFaceIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
#    self._integrators.append(nlfi)
    nlfi.thisown=0
%}
%pythonprepend BlockNonlinearForm::AddBdrFaceIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
    nlfi = args[0]
#    self._integrators.append(nlfi)
    nlfi.thisown=0
%}
}

LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementSpace *> &f,
    		        mfem::FiniteElementSpace *)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Array<int> *> &bdr_attr_is_ess,
 		        mfem::Array<int> *)
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::Vector *> &rhs, mfem::Vector *)

%include "fem/nonlinearform.hpp"

%extend mfem::NonlinearForm{
  mfem::SparseMatrix *GetGradientMatrix(const mfem::Vector &x) const
  {
      return  dynamic_cast<mfem::SparseMatrix *>(&$self->GetGradient(x));
  }
};

