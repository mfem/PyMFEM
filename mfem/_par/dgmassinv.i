//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") dgmassinv

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_DGMASSINV

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/exception.i"

%import "fespace.i"

%pythonprepend mfem::DGMassInverse::DGMassInverse %{
    # store fes and coeff passed as Python object
    self._link = args
%}

%include "fem/dgmassinv.hpp"

#endif //FILE_EXISTS_...
