%module(package="mfem._par") matrix

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>
#include "mfem.hpp"
#include "numpy/arrayobject.h"  
#include "../common/io_stream.hpp"  
#include "pyoperator.hpp"
%}

%begin %{
#define PY_SSIZE_T_CLEAN
%}

%init %{
import_array();
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
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(Matrix, Print)
#endif
