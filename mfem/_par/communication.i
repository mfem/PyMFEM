%module(package="mfem._par") communication
%feature("autodoc", "1");

%{
#include <fstream>
#include <iostream>  
#include <mpi.h>
#include "../common/io_stream.hpp"
#include "mfem.hpp"
#include "numpy/arrayobject.h"
%}

%include "../common/mfem_config.i"

%init %{
import_array();
%}

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
/*
%init %{
import_array();
%}
*/
%include "exception.i"
%import "array.i"
%import "table.i"
%import "sets.i"

%ignore mfem::MPI_Session;
%include "general/communication.hpp"

/*
   void Save(std::ostream &out) const;
   void PrintInfo(std::ostream &out = mfem::out) const;
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(GroupTopology, Save)
OSTREAM_ADD_DEFAULT_FILE(GroupCommunicator, PrintInfo)
#endif
