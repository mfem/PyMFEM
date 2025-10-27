//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module (package="mfem._par", director="1") intrules_cut

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "../common/exception.i"
%import "array.i"
%import "../common/numpy_int_typemap.i"
%import "mem_manager.i"

%import "coefficient.i"

%include "fem/intrules_cut.hpp"
