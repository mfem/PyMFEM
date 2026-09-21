//
// Copyright (c) 2020-2026, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") filteredsolver

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pysolvers.hpp"  
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_LINALG_FILTEREDSOLVER

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/exception.i"

%import "solvers.i"

%include "linalg/filteredsolver.hpp"

#endif //FILE_EXISTS_...
