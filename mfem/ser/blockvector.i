%module blockvector

%{
#include "linalg/blockvector.hpp"
#include "numpy/arrayobject.h"
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%include "linalg/blockvector.hpp"
