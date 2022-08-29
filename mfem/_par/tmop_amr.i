%module(package="mfem._par", directors="1") tmop_amr
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%include "exception.i"
%import "../common/exception_director.i"

%import tmop.i
%import estimators.i
%import nonlinearform.i
%import pnonlinearform.i

%include "fem/tmop_amr.hpp"
