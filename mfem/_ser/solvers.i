%module(package="mfem._ser") solvers

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
%ignore mfem::IterativeSolverMonitor::SetIterativeSolver;

%include "linalg/solvers.hpp"
