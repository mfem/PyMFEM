%module(package="mfem._ser") mesh_operators

%{
#include "numpy/arrayobject.h"
#include "mesh/mesh_operators.hpp"  
  %}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "mesh.i"
%import "estimators.i"
%import "../common/exception.i"

%include "mesh/mesh_operators.hpp"  
