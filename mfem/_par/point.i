%module(package="mfem._par") point

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "pyoperator.hpp"      
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%import "element.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

// to give index array as list
INT_TO_INTARRAY_IN(const int *ind)
INTARRAY_OUT_TO_INT(int *GetVertices)
INT_DEFAULT_NEGATIVE_ONE(int attr = -1)

%include "mesh/point.hpp"

