%module blockmatrix

%{
#include "linalg/blockmatrix.hpp"
#include "numpy/arrayobject.h"    
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "matrix.i"
%import "sparsemat.i"
%include "linalg/blockmatrix.hpp"
