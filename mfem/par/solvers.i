%module solvers
%{
#include <mpi.h>
#include "config/config.hpp"
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "config/_config.hpp" // include mfem MACRO
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "sparsemat.i"

%include "linalg/solvers.hpp"

