//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") vertex

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}
%include "exception.i"
%import "element.i"
%include "../common/exception.i"

%typemap(in) (char *str, int len) {
  $1 = PyString_AsString($input);
  $2 = PyString_Size($input);
};


%include "../common/deprecation.i"
DEPRECATED_OVERLOADED_METHOD(mfem::Vertex::SetCoords,
    	                     Vertex::SetCoords(const double *p) is deprecated,
			     len(args) == 2)
 //%feature("pythonprepend") mfem::Vertex::SetCoords%{
 //   if len(args) == 2:
 //       warnings.warn("SetCoords(const double *p) is deprecated",
 //   	              DeprecationWarning,)
 // %}
%include "mesh/vertex.hpp"

