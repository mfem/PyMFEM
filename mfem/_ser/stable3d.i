%module(package="mfem._ser") stable3d

%feature("autodoc", "1");

%{
#include <fstream>
#include  "mfem.hpp"
#include "../common/io_stream.hpp"
#include "pyoperator.hpp"        
#include "numpy/arrayobject.h"    
%}

%begin %{
#define PY_SSIZE_T_CLEAN
%}
%init %{
import_array();
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


