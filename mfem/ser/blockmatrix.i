%module blockmatrix

%{
#include "iostream_typemap.hpp"        
#include "linalg/blockmatrix.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"     
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "matrix.i"
%import "sparsemat.i"
%import "ostream_typemap.i"
%import "../common/ignore_common_functions.i"

%include "linalg/blockmatrix.hpp"
