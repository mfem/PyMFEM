//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") estimators_hdg

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pycoefficient.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_DARCY_ESTIMATORS_HDG

%init %{
import_array();
%}
%include "exception.i"
%include "../common/exception.i"

%import "bilininteg.i"

%include "fem/darcy/estimators_hdg.hpp"

#endif
