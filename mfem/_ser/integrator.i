//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") integrator

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"  
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_INTEGRATOR

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/exception.i"

%import "intrules.i"


%include "fem/integrator.hpp"

#endif
