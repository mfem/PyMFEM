%module geom
%{
#include "fem/geom.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
  
%import intrules.i
%import densemat.i

%immutable RefPts;
%immutable GlobGeometryRefiner;
%include "fem/geom.hpp"
