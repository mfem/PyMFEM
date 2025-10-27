//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser", directors="1") tmop_amr
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

%include "../common/mfem_config.i"

%include "exception.i"
%import "../common/exception_director.i"

%import tmop.i
%import nonlinearform.i
%import estimators.i

%include "fem/tmop_amr.hpp"
