//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") complexweakform
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/dpg/util/complexweakform.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pylininteg.hpp"  
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%init %{
import_array1(-1);
%}

%include "../common/existing_mfem_headers.i"

#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_COMPLEXWEAKFORM

%inline %{
#include "miniapps/dpg/util/complexweakform.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "gridfunc.i"
%import "mesh.i"
%import "solvers.i"
%import "operators.i"
%import "blockmatrix.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)

%include "miniapps/dpg/util/complexweakform.hpp"

#endif
