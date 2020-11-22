%module(package="mfem._ser") blockvector

%{
#include "linalg/blockvector.hpp"
#include "numpy/arrayobject.h"
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%include "exception.i"
%import "array.i"
%import "vector.i"
%import "../common/exception.i"

%import "mem_manager.i"

%include "linalg/blockvector.hpp"

