%module(package="mfem._ser") device
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

%include "../common/deprecation.i"
DEPRECATED_METHOD(mfem::Device::GetMemoryType())
DEPRECATED_METHOD(mfem::Device::GetMemoryClass())

%include "general/device.hpp"
