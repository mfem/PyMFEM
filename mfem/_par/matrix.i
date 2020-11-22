%module(package="mfem._par") matrix

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>
#include "numpy/arrayobject.h"  
#include "io_stream.hpp"  
#include "pyoperator.hpp"
#include "linalg/matrix.hpp"  
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
