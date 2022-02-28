/*

   constraints.i

*/
%module(package="mfem._par") constraints
%feature("autodoc", "1");
%{
#include "mfem.hpp"      
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pysolvers.hpp"
%}
%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "fespace.i"
%import "operators.i"
%import "sparsemat.i"
%import "hypre.i"
%import "solvers.i"

%include "linalg/constraints.hpp"
