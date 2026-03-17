//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") particlevector

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"  
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "vector.i"
%import "ordering.i"

%include "linalg/particlevector.hpp"

