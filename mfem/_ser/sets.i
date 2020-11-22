%module(package="mfem._ser") sets
%{
#include "general/sets.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%import "array.i"
%import "table.i"
%import "../common/exception.i"
%include "general/sets.hpp"
