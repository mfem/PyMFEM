//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") sets
%{
#include "general/sets.hpp"
#include "numpy/arrayobject.h"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "array.i"
%import "table.i"
%import "../common/exception.i"
%include "general/sets.hpp"
