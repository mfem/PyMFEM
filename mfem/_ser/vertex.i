%module(package="mfem._ser") vertex
%{
#include  "mfem.hpp"
#include "mesh/vertex.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%import "element.i"
%include "../common/exception.i"

%include "mesh/vertex.hpp"

