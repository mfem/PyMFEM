//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") fe_nurbs
%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}
%include "exception.i"
%import "fe_base.i"
%import "element.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"


%include "fem/fe/fe_nurbs.hpp"
