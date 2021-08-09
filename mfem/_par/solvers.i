%module(package="mfem._par", directors="1")  solvers
%{
#include <mpi.h>
#include "config/config.hpp"
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"
using namespace mfem;
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "sparsemat.i"
%import "../common/exception.i"

%ignore mfem::IterativeSolverMonitor::SetIterativeSolver;
%feature("director") mfem::IterativeSolverMonitor;
%include "linalg/solvers.hpp"
