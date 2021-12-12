%module(package="mfem._ser") mem_manager
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
%include "../common/enum_class.i"

//%include "general/mem_manager.hpp"
%include "../../headers/mem_manager.hpp"

ENUM_CLASS_WRAP(_mem_manager, MemoryType)
ENUM_CLASS_WRAP(_mem_manager, MemoryClass)
