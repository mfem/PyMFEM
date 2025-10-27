//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par", directors="1") nonlininteg
%{
#include "mfem.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pynonlininteg.hpp"
using namespace mfem;
%}
/*
%init %{
import_array1(-1);
%}
*/

%include "exception.i"
%import vector.i
%import operators.i
%import fespace.i
%import eltrans.i
%import integrator.i
%import "../common/exception_director.i"

%include "fem/nonlininteg.hpp"

%feature("director") mfem::PyNonlinearFormIntegrator;
%include "../common/pynonlininteg.hpp"

