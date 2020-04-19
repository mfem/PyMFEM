%module(package="mfem._ser") vtk
%{
#include  "mfem.hpp"
#include "mesh/vtk.hpp"
#include "numpy/arrayobject.h"    
#include "io_stream.hpp"
%}

%init %{
import_array();
%}
%include "exception.i"
%import "element.i"
%include "../common/exception.i"

%import "../common/numpy_int_typemap.i"
%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%include "mesh/vtk.hpp"

