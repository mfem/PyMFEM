//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") ode
%{
#include  "mfem.hpp"
#include "linalg/ode.hpp"
#include "../common/pyoperator.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
%}

%init %{
import_array1(-1);
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

