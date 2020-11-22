%module(package="mfem._par") globals
%{
#include  "mfem.hpp"
#include "general/globals.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%include "std_string.i"

%immutable out;
%immutable err;

%include "general/globals.hpp"
