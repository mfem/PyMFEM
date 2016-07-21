%module plinearform
%{
#include <mpi.h>
#include "fem/plinearform.hpp"
#include "numpy/arrayobject.h"
#define MFEM_USE_MPI  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}

%import linearform.i
%import pfespace.i
%import hypre.i
#define MFEM_USE_MPI
%include "fem/plinearform.hpp"
