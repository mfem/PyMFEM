%module(package="mfem._par") prestriction
%{
#include  "mfem.hpp"
#include  "fem/prestriction.hpp"
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

%import "restriction.i"

%include "fem/prestriction.hpp"

