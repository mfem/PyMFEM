%module(package="mfem._par") mem_manager
%{
#include  "mfem.hpp"
#include "general/mem_manager.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%include "../../headers/mem_manager.hpp"

