%module vertex
%{
#include "mesh/vertex.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%include "mesh/vertex.hpp"

