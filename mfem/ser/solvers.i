%module solvers

%{
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
%}

%import "linalg/vector.hpp"
%import "operators.i"
%import "linalg/matrix.hpp"
%import "linalg/sparsemat.hpp"
%include "linalg/solvers.hpp"
