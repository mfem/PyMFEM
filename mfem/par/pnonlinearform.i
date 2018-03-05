%module pnonlinearform
%{
#include <mpi.h>
#include "config/config.hpp"
#include "fem/pnonlinearform.hpp"
#include "fem/linearform.hpp"
#include "numpy/arrayobject.h"  
#include "pyoperator.hpp"           
%}

%include "config/config.hpp" // include mfem MACRO
#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

/*
%init %{
import_array();
%}
*/
%include "exception.i"
%import "vector.i"
%import "nonlinearform.i"
%import "blockoperator.i"
%import "pfespace.i"
%import "pgridfunc.i"
%import "../common/exception.i"

%pointer_class(int, intp);

%include "fem/pnonlinearform.hpp"
