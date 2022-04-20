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

//  enforce that device class is created only once
%pythoncode %{
  singleton_device = None
%}
%feature("shadow") mfem::Device::__new__ %{
    def __new__(cls, *args, **kwargs):
        if globals()["singleton_device"] is None:  
             instance =  super(Device, cls).__new__(cls)
             globals()["singleton_device"] = instance
        return globals()["singleton_device"]
%}
%pythonappend mfem::Device::Device %{
  globals()["singleton_device"] = self
%}

%include "general/device.hpp"

%extend mfem::Device{
  // placeholder to create __new__  
  void __new__(void){}
};

