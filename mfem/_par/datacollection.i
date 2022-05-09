%module(package="mfem._par") datacollection
%{
#include "mfem.hpp"
#include "pyoperator.hpp"    
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"      
%}

%include "../common/mfem_config.i"

%init %{
import_array();
%}

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "globals.i"
%import "mesh.i"
%import "gridfunc.i"
%import "pgridfunc.i"

%include "fem/datacollection.hpp"
