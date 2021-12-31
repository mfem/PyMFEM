%module(package="mfem._par") segment
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

LIST_TO_INTARRAY_IN(const int *ind, 2)
INTARRAY_OUT_TO_TUPLE(int *GetVertices, 2)

%include "mesh/segment.hpp"

