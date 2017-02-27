%module ncmesh
%{
#include "iostream_typemap.hpp"           
#include "mesh/ncmesh.hpp"
#include "numpy/arrayobject.h"    
%}

%import "mesh.i"
%import "array.i"
%import "fem/geom.hpp"

%import "ostream_typemap.i"

%immutable embeddings;

%include  "mesh/ncmesh.hpp"
