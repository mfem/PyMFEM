%module sparsesmoothers
%{
#include "linalg/sparsesmoothers.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%import "vector.i"
%import "operators.i"
%import "sparsemat.i"
%import "matrix.i"

%include "linalg/sparsesmoothers.hpp" 
