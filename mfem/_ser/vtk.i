%module(package="mfem._ser") vtk
%{
#include  "mfem.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
#include "../common/io_stream.hpp"
%}
%begin %{
#define PY_SSIZE_T_CLEAN
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
