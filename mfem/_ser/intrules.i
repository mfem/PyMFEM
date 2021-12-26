%module (package="mfem._ser") intrules

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

%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::IntegrationPoint)

%import "../common/array_instantiation_macro.i"
IGNORE_OBJ_METHODS(IntegrationPoint)
INSTANTIATE_ARRAY(IntegrationPoint)

%include "fem/intrules.hpp"
