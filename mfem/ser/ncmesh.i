%module ncmesh
%{
#include "iostream_typemap.hpp"           
#include "mesh/ncmesh.hpp"
#include "numpy/arrayobject.h"    
%}

%include "exception.i"
%include "../common/cpointers.i"
%import "mesh.i"
%import "array.i"
%import "fem/geom.hpp"
%import "ostream_typemap.i"
%import "../common/exception.i"

%immutable embeddings;

%include  "mesh/ncmesh.hpp"
