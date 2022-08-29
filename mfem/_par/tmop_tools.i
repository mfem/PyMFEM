%module(package="mfem._par", directors="1") tmop_tools
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pysolvers.hpp"  
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

%include "../common/typemap_macros.i"

%import tmop.i
%import bilinearform.i
%import pbilinearform.i
%import solvers.i

%include "fem/tmop_tools.hpp"
