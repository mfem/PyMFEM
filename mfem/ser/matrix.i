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

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "array.i"
%import "ostream_typemap.i"
%import "../common/exception.i"

%include "linalg/matrix.hpp"
