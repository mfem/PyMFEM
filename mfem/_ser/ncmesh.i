%module(package="mfem._ser") ncmesh

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>
#include "io_stream.hpp"           
#include "mesh/ncmesh.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
   import_array();
%}

%include "exception.i"
%include "../common/cpointers.i"
%import "mesh.i"
%import "array.i"
%import "fem/geom.hpp"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%immutable embeddings;

%include  "mesh/ncmesh.hpp"

#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_STDOUT_FILE(NCMesh, PrintVertexParents)
OSTREAM_ADD_DEFAULT_STDOUT_FILE(NCMesh, PrintCoarseElements)
OSTREAM_ADD_DEFAULT_FILE(NCMesh, PrintStats)
#endif
