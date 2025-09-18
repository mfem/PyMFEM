//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") hash

%{
#include "general/hash.hpp"
#include "numpy/arrayobject.h"
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%include "exception.i"
%import "array.i"
%import "vector.i"
%import "../common/exception.i"


%include "general/hash.hpp"

