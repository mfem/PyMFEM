%module communication
%{
#include <mpi.h>
#include "iostream_typemap.hpp"      
#include "config/config.hpp"    
#include "general/sets.hpp"
#include "general/communication.hpp"
#include "numpy/arrayobject.h"
%}

%include "../common/mfem_config.i"

%init %{
import_array();
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

%import ostream_typemap.i

%include "general/communication.hpp"
