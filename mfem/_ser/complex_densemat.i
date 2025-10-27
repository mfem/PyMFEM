//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") complex_densemat
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "array.i"
%import "vector.i"
%import "densemat.i"
%import "complex_operator.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

%include "linalg/complex_densemat.hpp"
