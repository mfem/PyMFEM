%module(package="mfem._par") device
%{
#include  "mfem.hpp"
#include "general/device.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "std_string.i"

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"



%import "mem_manager.i"

%include "general/device.hpp"
