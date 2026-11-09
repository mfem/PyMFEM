//
// Copyright (c) 2020-2026, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par")  sidredatacollection
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_SIDREDATACOLLECTION

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "datacollection.i"

%include "fem/sidredatacollection.hpp"

#endif
