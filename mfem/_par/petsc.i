%module(package="mfem._par") petsc
%{
#include "mfem.hpp"

#include "../common/pycoefficient.hpp"  
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"  
#include "../common/io_stream.hpp"  
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
%import "fe.i"
%import "fe_fixed_order.i"
%import "element.i"
%import "mesh.i"
%import "operators.i"
%import "coefficient.i"

%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%import   "petscconf.h"
%include "linalg/petsc.hpp"





