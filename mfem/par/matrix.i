%module matrix

%{
#include "linalg/matrix.hpp"
#include "pyoperator.hpp"
#include "iostream_typemap.hpp"  
%}

%import "vector.i"
%import "operators.i"
%import "array.i"
%import "ostream_typemap.i"

%include "linalg/matrix.hpp"
