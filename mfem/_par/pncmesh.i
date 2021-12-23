%module(package="mfem._par") pncmesh

%feature("autodoc", "1");

%{
#include <mpi.h>
#include "../common/io_stream.hpp"       
#include "config/config.hpp"
#include "mesh/mesh_headers.hpp"
#include "mpi4py/mpi4py.h"
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
%import "mesh.i"
 //%import "ncmesh.i"
%import "communication.i"

%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%pointer_class(int, intp);

%include "mesh/pncmesh.hpp"

/*
  void Dump(std::ostream &os) const;
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_STDOUT_FILE(ElementSet, Dump)
#endif
