%module matrix

%{
#include "linalg/matrix.hpp"
#include "iostream_typemap.hpp"      
#include "pyoperator.hpp"             
%}

%import "vector.i"
%import "operators.i"
%import "array.i"
%import "ostream_typemap.i"

%include "linalg/matrix.hpp"
