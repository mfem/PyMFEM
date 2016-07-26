%module solvers

%{
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
#include "pyoperator.hpp"     
%}

%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "sparsemat.i"

%include "linalg/solvers.hpp"
