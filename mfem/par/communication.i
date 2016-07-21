%module communication
%{
#include <mpi.h>
#include "general/sets.hpp"
#include "general/communication.hpp"
#define MFEM_USE_MPI  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
/*
%init %{
import_array();
%}
*/
%import array.i
%import table.i
%import sets.i
#define MFEM_USE_MPI  
%include "general/communication.hpp"
