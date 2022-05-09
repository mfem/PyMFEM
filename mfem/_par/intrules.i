%module (package="mfem._par") intrules

%{
#include "fem/intrules.hpp"
#include "numpy/arrayobject.h"
%}

%init %{
import_array();
%}

%include "exception.i"
%import "../common/exception.i"
%import "array.i"
%import "../common/numpy_int_typemap.i"
%import "mem_manager.i"

%immutable IntRules;
%immutable RefinedIntRules;

/* define IntegrationPointArray and IntegrationRulePtrArray */
%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::IntegrationPoint, 0)
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::IntegrationRule *, 1)

%import "../common/data_size_typemap.i"
XXXPTR_SIZE_IN(mfem::IntegrationPoint *data_, int asize, mfem::IntegrationPoint)
XXXPTR_SIZE_IN(mfem::IntegrationRule **data_, int asize, mfem::IntegrationRule *)

%import "../common/array_instantiation_macro.i"
IGNORE_ARRAY_METHODS(mfem::IntegrationPoint)
INSTANTIATE_ARRAY(IntegrationPoint)
IGNORE_ARRAY_METHODS(mfem::IntegrationRule *)
INSTANTIATE_ARRAY0(IntegrationRule *, IntegrationRule, 1)

%include "fem/intrules.hpp"

