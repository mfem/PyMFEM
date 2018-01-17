%module ode
%{
#include "linalg/ode.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%import "vector.i"
%import "array.i"
%import "operators.i"

%typemap(in) double &t (double temp){
  temp3 = PyFloat_AsDouble($input);
  $1 = &temp3;
 }
%typemap(in) double &dt (double dtemp){
  dtemp4 = PyFloat_AsDouble($input);
  $1 = &dtemp4;
}
%typemap(argout) double &t {
  %append_output(PyFloat_FromDouble(*$1));
}
%typemap(argout) double &dt {
  %append_output(PyFloat_FromDouble(*$1));
 }


%include "linalg/ode.hpp"

