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

%ignore mfem::Array<mfem::IntegrationPoint>::Union;
%ignore mfem::Array<mfem::IntegrationPoint>::Find;
%ignore mfem::Array<mfem::IntegrationPoint>::FindSorted;
%ignore mfem::Array<mfem::IntegrationPoint>::Sort;
%ignore mfem::Array<mfem::IntegrationPoint>::DeleteFirst;
%ignore mfem::Array<mfem::IntegrationPoint>::Unique;
%ignore mfem::Array<mfem::IntegrationPoint>::PartialSum;
%ignore mfem::Array<mfem::IntegrationPoint>::Sum;
%ignore mfem::Array<mfem::IntegrationPoint>::IsSorted;
%ignore mfem::Array<mfem::IntegrationPoint>::Save;
%ignore mfem::Array<mfem::IntegrationPoint>::Max;
%ignore mfem::Array<mfem::IntegrationPoint>::Min;
%ignore mfem::Array<mfem::IntegrationPoint>::Print;
%ignore mfem::Array<mfem::IntegrationPoint>::Load;
%template(IntegrationPointArray) mfem::Array<mfem::IntegrationPoint>;

%include "fem/intrules.hpp"

