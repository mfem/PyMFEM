//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") particleset

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "../common/io_stream.hpp"  
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"
%include "../common/kernel_dispatch.i"

%import "array.i"
%import "ordering.i"
%import "particlevector.i"

%include "fem/particleset.hpp"

