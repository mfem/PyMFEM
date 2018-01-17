%module communication
%{
#include <mpi.h>
#include "iostream_typemap.hpp"      
#include "config/config.hpp"    
#include "general/sets.hpp"
#include "general/communication.hpp"
#include "numpy/arrayobject.h"
%}

%include  "config/_config.hpp" // include mfem MACRO
%init %{
import_array();
%}

%include "config/_config.hpp" // include mfem MACRO
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
