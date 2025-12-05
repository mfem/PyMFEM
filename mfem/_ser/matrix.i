//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") matrix

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
%}

%begin %{
#define PY_SSIZE_T_CLEAN
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "array.i"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%include "linalg/matrix.hpp"

/*
  virtual void Print (std::ostream & out = mfem::out, int width_ = 4) const;
*/
OSTREAM_ADD_DEFAULT_FILE(Matrix, Print)
 //OSTREAM_ADD_DEFAULT_FILE_ARG1(Matrix, Print, int width_=4, width_)
