//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") mesh_operators

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
  %}
// initialization required to return numpy array from SWIG
%init %{
import_array1(-1);
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "mesh.i"
%import "estimators.i"
%import "../common/exception.i"

%include "../common/typemap_macros.i"
LIST_TO_MFEMOBJ_POINTERARRAY_IN(mfem::IntegrationRule const *irs_[],  mfem::IntegrationRule *, 1)

%include "mesh/mesh_operators.hpp"
