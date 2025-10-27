//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") stable3d

%feature("autodoc", "1");

%{
#include <fstream>
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
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%include "general/stable3d.hpp"

/*
void Print(std::ostream &out = mfem::out) const;
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(STable3D, Print)
#endif


