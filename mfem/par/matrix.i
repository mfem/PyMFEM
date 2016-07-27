%module matrix

%{
#include "linalg/matrix.hpp"
#include "pyoperator.hpp"               
%}
%import "vector.i"
%import "operators.i"
%import "array.i"
%include "linalg/matrix.hpp"
