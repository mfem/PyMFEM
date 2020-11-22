%module(package="mfem._par") sparsesmoothers
%{
#include "linalg/sparsesmoothers.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "sparsemat.i"
%import "matrix.i"
%import "../common/exception.i"

%include "linalg/sparsesmoothers.hpp" 
