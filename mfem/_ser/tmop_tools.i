//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser", directors="1") tmop_tools
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"

%include "exception.i"
%import "../common/exception_director.i"

%include "../common/typemap_macros.i"

#ifdef MFEM_USE_GSLIB
%import gslib.i
#endif
%import tmop.i
%import bilinearform.i
%import solvers.i

%include "fem/tmop_tools.hpp"
