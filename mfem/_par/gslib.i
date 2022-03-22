%module(package="mfem._par") gslib
%{
#include "mfem.hpp"
#include "fem/gslib.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"  
#include "../common/pycoefficient.hpp"
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array();
%}

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"
%import vector.i
%import mesh.i
%import gridfunc.i

%include "fem/gslib.hpp"
