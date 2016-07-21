%module fe

%{
#include "config/config.hpp"
#include "linalg/linalg.hpp"
#include "fem/intrules.hpp"
#include "fem/fe.hpp"
%}

%immutable;
%ignore poly1d;
%mutable;

%import "array.i"
%import "vector.i"
%import "intrules.i"
%import "densemat.i"
%import "sparsemat.i"

%include "fem/fe.hpp"


