%module(package="mfem._ser") ncmesh

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>
#include "../common/io_stream.hpp"
#include "general/zstr.hpp"
#include "general/globals.hpp"
#include "mesh/ncmesh.hpp"
#include "fem/gridfunc.hpp"
#include "fem/bilinearform.hpp"    
#include "numpy/arrayobject.h"
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
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::Refinement)

%immutable embeddings;

%include  "mesh/ncmesh.hpp"

#ifndef SWIGIMPORTED
 //OSTREAM_ADD_DEFAULT_STDOUT_FILE(NCMesh, PrintVertexParents)
 //OSTREAM_ADD_DEFAULT_STDOUT_FILE(NCMesh, PrintCoarseElements)
 //OSTREAM_ADD_DEFAULT_FILE(NCMesh, PrintStats)
#endif

namespace mfem{
   %ignore Array<Refinement>::Union;
   %ignore Array<Refinement>::Find;
   %ignore Array<Refinement>::FindSorted;
   %ignore Array<Refinement>::Sort;
   %ignore Array<Refinement>::DeleteFirst;
   %ignore Array<Refinement>::Unique;
   %ignore Array<Refinement>::PartialSum;
   %ignore Array<Refinement>::Sum;
   %ignore Array<Refinement>::IsSorted;
   %ignore Array<Refinement>::Save;
   %ignore Array<Refinement>::Max;
   %ignore Array<Refinement>::Min;
   %ignore Array<Refinement>::Print;
   %ignore Array<Refinement>::Load;
}
%import "../common/array_instantiation_macro.i"
INSTANTIATE_ARRAY(Refinement)

