%module pnonlinearform
%{
#include <mpi.h>
#include "fem/pnonlinearform.hpp"
#include "fem/linearform.hpp"  
#define MFEM_USE_MPI  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
/*
%init %{
import_array();
%}
*/
%import nonlinearform.i
%import pfespace.i
%import pgridfunc.i
#define MFEM_USE_MPI
%include "fem/pnonlinearform.hpp"
