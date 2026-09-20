//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") lor

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_LOR_LOR

%init %{
import_array1(-1);
%}
%include "exception.i"
%include "../common/exception.i"

%import "pbilinearform.i"
%import "pfespace.i"

%include "fem/lor/lor.hpp"

%template(LORHypreBoomerAMG) mfem::LORSolver<mfem::HypreBoomerAMG>;
%template(LORHypreAMS) mfem::LORSolver<mfem::HypreAMS>;
%template(LORHypreADS) mfem::LORSolver<mfem::HypreADS>;
#endif //FILE_EXISTS_...
