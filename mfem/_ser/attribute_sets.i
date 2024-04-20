%module(package="mfem._ser") attribute_sets
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MESH_ATTRIBUTE_SETS

%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%import "array.i"
%import "arrays_by_name.i"
%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%include "mesh/attribute_sets.hpp"

#endif
