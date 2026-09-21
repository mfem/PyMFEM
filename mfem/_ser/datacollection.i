//
// Copyright (c) 2020-2026, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") datacollection
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "globals.i"
%import "mesh.i"
%import "gridfunc.i"

%ignore mfem::ParaViewDataCollection::GetCoeffFieldMap;
%ignore mfem::ParaViewDataCollection::GetVCoeffFieldMap;

%include "fem/datacollection.hpp"
