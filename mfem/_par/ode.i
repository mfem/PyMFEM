%module(package="mfem._par") ode
%{
#include  "mfem.hpp"
#include "linalg/ode.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
#include "../common/io_stream.hpp"
%}

%init %{
import_array();
%}
%include "exception.i"
%import "vector.i"
%import "array.i"
%import "operators.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%typemap(in) double &t (double temp){
  temp = PyFloat_AsDouble($input);
  $1 = &temp;
 }
%typemap(in) double &dt (double dtemp){
  dtemp = PyFloat_AsDouble($input);
  $1 = &dtemp;
}
%typemap(argout) double &t {
  %append_output(PyFloat_FromDouble(*$1));
}
%typemap(argout) double &dt {
  %append_output(PyFloat_FromDouble(*$1));
 }


%include "linalg/ode.hpp"

