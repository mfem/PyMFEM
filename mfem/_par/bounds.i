//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") bounds

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"  
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_BOUNDS

%init %{
import_array1(-1);
%}
%include "exception.i"
%import "vector.i"
%import "fespace.i"
%include "../common/exception.i"

%include "fem/bounds.hpp"

#endif
