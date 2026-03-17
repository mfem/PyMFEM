//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") config

%{
#include "mfem.hpp"
%}


%include "../common/mfem_config.i"

#ifdef MFEM_USE_SINGLE
%pythoncode %{
  MFEM_USE_SINGLE = 1
%}
#elif defined MFEM_USE_DOUBLE
%pythoncode %{
  MFEM_USE_SINGLE = 0
%}
#else
#error "Either DOUBLE or SINGLE precision must be specified"
#endif


