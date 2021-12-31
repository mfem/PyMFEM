%module(package="mfem._ser") point

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

 //
%include "../common/deprecation.i"
DEPRECATED_OVERLOADED_METHOD(mfem::Point::GetNFaces,
    	                     Point::GetNFaces(int & nFaceVertices) is deprecated,
			     len(args) == 1)

%include "mesh/point.hpp"

