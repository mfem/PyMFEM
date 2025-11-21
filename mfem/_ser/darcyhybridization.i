//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") darcyhybridization

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pycoefficient.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_DARCY_DARCYHYBRIDIZATION

%init %{
import_array();
%}
%include "exception.i"
%include "../common/exception.i"

 //%ignore mfem::DarcyHybridization::m_nlfi_u;
%ignore mfem::DarcyHybridization;
%import "nonlininteg.i"

%include "fem/darcy/darcyhybridization.hpp"

#endif
