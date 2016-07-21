%module solvers
%{
#include <mpi.h>  
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
#define MFEM_USE_MPI    
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "sparsemat.i"
#define MFEM_USE_MPI  
%include "linalg/solvers.hpp"

