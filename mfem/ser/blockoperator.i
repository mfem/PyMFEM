%module blockoperator

%{
#include "linalg/blockoperator.hpp"
#include "numpy/arrayobject.h"    
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "operators.i"
%include "linalg/blockoperator.hpp"
