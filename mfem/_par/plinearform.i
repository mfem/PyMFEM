//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") plinearform
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "linearform.i"
%import "pfespace.i"
%import "pgridfunc.i"
%import "hypre.i"
%import "../common/exception.i"

%newobject mfem::ParLinearForm::ParallelAssemble;

%include "fem/plinearform.hpp"

