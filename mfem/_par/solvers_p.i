%module(package="mfem._par") solvers_p

%{
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
%}
#define MFEM_USE_MPI  1
%import "linalg/vector.hpp"
%import "operators.i"
%import "linalg/matrix.hpp"
%import "linalg/sparsemat.hpp"
%include "linalg/solvers.hpp"
