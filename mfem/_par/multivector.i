//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") multivector

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_LINALG_MULTIVECTOR

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/exception.i"

%import "array.i"
%import "vector.i"

%include "linalg/multivector.hpp"

#endif //FILE_EXISTS_...
