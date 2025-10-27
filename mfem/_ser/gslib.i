//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") gslib
%{
#include "mfem.hpp"
#include "fem/gslib.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/mfem_config.i"

%init %{
import_array1(-1);
%}

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"
%import array.i
%import vector.i
%import mesh.i
%import gridfunc.i

%include "fem/gslib.hpp"
