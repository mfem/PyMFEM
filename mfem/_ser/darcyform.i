//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") darcyform

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_DARCY_DARCYFORM

%init %{
import_array();
%}
%include "exception.i"
%include "../common/exception.i"

%import "fespace.i"
%import "blockvector.i"
%import "bilinearform.i"
%import "nonlinearform.i"

%include "fem/darcy/darcyform.hpp"

#endif

