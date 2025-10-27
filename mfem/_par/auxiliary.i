//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module (package="mfem._par") auxiliary

%{
#include <iostream>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "../common/mfem_config.i"
%include "exception.i"
%import "../common/exception.i"

%import "solvers.i"
%import "pfespace.i"

%include "linalg/auxiliary.hpp"


