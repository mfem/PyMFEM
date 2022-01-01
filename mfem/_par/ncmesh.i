%module(package="mfem._par") ncmesh

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>
#include "mfem.hpp"
#include "../common/io_stream.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"  
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}

%include "exception.i"
%import "globals.i"
%import "mesh.i"
%import "array.i"
%import "fem/geom.hpp"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::Refinement, 0)

%import "../common/data_size_typemap.i"
XXXPTR_SIZE_IN(mfem::Refinement *data_, int asize, mfem::Refinement)

%immutable embeddings;
%include  "mesh/ncmesh.hpp"

 /*
  void PrintVertexParents(std::ostream &out) const;
  void PrintCoarseElements(std::ostream &out) const;
  void PrintStats(std::ostream &out = mfem::out) const;
 */

#ifndef SWIGIMPORTED
 //OSTREAM_ADD_DEFAULT_STDOUT_FILE(NCMesh, PrintVertexParents)
 //OSTREAM_ADD_DEFAULT_STDOUT_FILE(NCMesh, PrintCoarseElements)
 //OSTREAM_ADD_DEFAULT_FILE(NCMesh, PrintStats)
#endif

%import "../common/array_instantiation_macro.i"
IGNORE_ARRAY_METHODS(mfem::Refinement)
INSTANTIATE_ARRAY(Refinement)

