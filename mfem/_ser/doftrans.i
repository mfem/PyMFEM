//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") doftrans
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
%import "vector.i"
%import "densemat.i"
%import "intrules.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"


%include "fem/doftrans.hpp"
