%module vertex
%{
#include "mesh/vertex.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%import "element.i"
%include "mesh/vertex.hpp"

