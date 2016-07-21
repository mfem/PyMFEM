%module ncmesh
%{
  #include "mesh/ncmesh.hpp"
  #include "numpy/arrayobject.h"    
%}

%import "mesh.i"
%import "array.i"
%import "fem/geom.hpp"

%immutable embeddings;

%include  "mesh/ncmesh.hpp"
