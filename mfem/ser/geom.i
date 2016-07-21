%module geom
%{
#include "fem/geom.hpp"
%}
%import intrules.i
%import densemat.i

%immutable RefPts;
%include "fem/geom.hpp"
