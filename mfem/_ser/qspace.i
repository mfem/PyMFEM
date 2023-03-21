%module(package="mfem._ser") qspace
%{
#include "mfem.hpp"
#include "../common/io_stream.hpp"  
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
#include "numpy/arrayobject.h"    
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_QSPACE

%init %{
import_array();
%}
%include "exception.i"
%import "fe.i"
%import "fe_fixed_order.i"
%import "element.i"
%import "mesh.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%include "fem/qspace.hpp"

OSTREAM_ADD_DEFAULT_STDOUT_FILE(QuadratureSpace, Save)

#endif



