%module(package="mfem._par") ncmesh

%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>

#include "../common/io_stream.hpp"
#include "mesh/ncmesh.hpp"
#include "general/zstr.hpp"
#include "fem/gridfunc.hpp"
#include "fem/bilinearform.hpp"  
#include "numpy/arrayobject.h"      
%}

%init %{
import_array();
%}

%include "exception.i"
%import "mesh.i"
%import "array.i"
%import "fem/geom.hpp"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

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

namespace mfem{
  %ignore Array<Refinement>::Union;  //
   %ignore Array<Refinement>::Find;
   %ignore Array<Refinement>::FindSorted;
   %ignore Array<Refinement>::Sort; //
   %ignore Array<Refinement>::DeleteFirst; //
   %ignore Array<Refinement>::Unique;  //
   %ignore Array<Refinement>::PartialSum;
   %ignore Array<Refinement>::Sum;
   %ignore Array<Refinement>::IsSorted;
   %ignore Array<Refinement>::Save;
   %ignore Array<Refinement>::Max;
   %ignore Array<Refinement>::Min;
   %ignore Array<Refinement>::Print;
   %ignore Array<Refinement>::Load;

   %template(RefinementArray) Array<Refinement>;
}
