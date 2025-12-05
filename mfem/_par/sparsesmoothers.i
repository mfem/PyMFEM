//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") sparsesmoothers
%{
#include "linalg/sparsesmoothers.hpp"
#include "../common/pyoperator.hpp"
#include "numpy/arrayobject.h"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "sparsemat.i"
%import "matrix.i"
%import "../common/exception.i"

%include "linalg/sparsesmoothers.hpp"
