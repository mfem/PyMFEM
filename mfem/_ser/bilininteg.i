//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser", directors="1")  bilininteg
%{
#include "mfem.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/pylininteg.hpp"
#include "numpy/arrayobject.h"
  //using namespace mfem;
%}

%init %{
import_array1(-1);
%}

%include "exception.i"

%import "globals.i"
%import "array.i"
%import "coefficient.i"
%import "matrix.i"
%import "vector.i"
%import "gridfunc.i"
%import "fespace.i"
%import "fe_coll.i"
%import "intrules.i"
%import "densemat.i"
%import "sparsemat.i"
%import "lininteg.i"
%import "eltrans.i"
%import "linearform.i"
%import "fe.i"
%import "nonlininteg.i"
%include "../common/exception_director.i"

%include "../common/bilininteg_ext.i"

%ignore  mfem::MassIntegrator::SetupPA;

%include "../common/kernel_dispatch.i"
%include "fem/bilininteg.hpp"

%feature("director") mfem::PyBilinearFormIntegrator;
%include "../common/pybilininteg.hpp"
