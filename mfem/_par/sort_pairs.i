//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") sort_pairs
%{
#include  "mfem.hpp"
#include "general/globals.hpp"
#include "numpy/arrayobject.h"
%}

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%include "general/sort_pairs.hpp"

namespace mfem{
%template(intintintTriple) Triple<int, int, int>;
 }
