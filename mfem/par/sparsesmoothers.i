%module sparsesmoothers
%{
#include "linalg/sparsesmoothers.hpp"
%}

%import "vector.i"
%import "operators.i"
%import "sparsemat.i"
%import "matrix.i"

%include "linalg/sparsesmoothers.hpp" 
