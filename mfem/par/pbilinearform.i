%module pbilinearform
%{
#include <mpi.h>
#include "fem/pbilinearform.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"           
#define MFEM_USE_MPI  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);


%init %{
import_array();
%}

%import bilinearform.i
%import pfespace.i
%import hypre.i

#define MFEM_USE_MPI  
%include "fem/pbilinearform.hpp"
