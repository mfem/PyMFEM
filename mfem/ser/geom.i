%module geom
%{
#include "fem/geom.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%import "intrules.i"
%import "densemat.i"
%import "../common/exception.i"

%immutable RefPts;
%immutable GlobGeometryRefiner;
%include "fem/geom.hpp"
