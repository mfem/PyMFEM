//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") prestriction
%{
#include  "mfem.hpp"
#include "fem/prestriction.hpp" // note:prestrition.hpp is not inclured from mfem.hpp
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}
%include "exception.i"
%import "element.i"
%include "../common/exception.i"
%import "../common/numpy_int_typemap.i"

%import "restriction.i"

%include "fem/prestriction.hpp"

