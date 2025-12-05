//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") vertex

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%include "mesh/vertex.hpp"

