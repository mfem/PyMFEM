%module(package="mfem._ser") restriction
%{
#include  "mfem.hpp"
#include  "fem/restriction.hpp"
#include "numpy/arrayobject.h"  
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}
%include "exception.i"
%import "element.i"
%include "../common/exception.i"
%import "../common/numpy_int_typemap.i"

%import "operators.i"
%import "mesh.i"

%include "fem/restriction.hpp"

