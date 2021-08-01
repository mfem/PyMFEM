%module(package="mfem._ser", directors="1") solvers

%{
#include "linalg/handle.hpp"  
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "sparsemat.i"
%import "../common/exception_director.i"
%import "../common/operator_ptr_typemap.i"
%import "../common/exception_director.i"
%ignore mfem::IterativeSolverMonitor::SetIterativeSolver;
%feature("director") mfem::IterativeSolverMonitor;
%include "linalg/solvers.hpp"

