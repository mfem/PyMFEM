//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
/*

   constraints.i

*/
%module(package="mfem._ser") constraints
%feature("autodoc", "1");
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pyintrules.hpp"
%}
%init %{
import_array1(-1);
%}

%include "exception.i"
%import "vector.i"
%import "fespace.i"
%import "operators.i"
%import "sparsemat.i"
%import "solvers.i"

%include "linalg/constraints.hpp"
