//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") point

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}
%include "exception.i"
%import "element.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

// to give index array as list
INT_TO_INTARRAY_IN(const int *ind)
INTARRAY_OUT_TO_INT(int *GetVertices)
INT_DEFAULT_NEGATIVE_ONE(int attr = -1)

%include "mesh/point.hpp"

