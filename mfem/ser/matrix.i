%module matrix

%{
#include "linalg/matrix.hpp"
#include "iostream_typemap.hpp"      
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%import "vector.i"
%import "operators.i"
%import "array.i"
%import "ostream_typemap.i"

%include "linalg/matrix.hpp"
