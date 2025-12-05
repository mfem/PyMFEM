//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") hybridization
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
%}

%init %{
import_array1(-1);
%}

%import "handle.i"
%import "fespace.i"
%import "bilininteg.i"
%import "hypre.i"

%include "fem/hybridization.hpp"
