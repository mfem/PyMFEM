//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") vtk
%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}
%begin %{
#define PY_SSIZE_T_CLEAN
%}
%init %{
import_array1(-1);
%}
%include "exception.i"
%import "element.i"
%include "../common/exception.i"

%import "../common/numpy_int_typemap.i"
%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%include "mesh/vtk.hpp"
