//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
// deprecate an overload method
%define DEPRECATED_OVERLOADED_METHOD(method, message, condition)
%feature("pythonprepend") method %{
  if condition:
       import warnings
       warnings.warn("message is deprecated",
   	              DeprecationWarning,)
%}
%enddef
// deprecate an overload method
%define DEPRECATED_METHOD(method)
%feature("pythonprepend") method %{
   import warnings
   warnings.warn("method is deprecated",
                 DeprecationWarning,)
%}
%enddef  
