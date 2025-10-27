//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") pncmesh

%feature("autodoc", "1");

%{
#include <mpi.h>
#include "mpi4py/mpi4py.h"
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "mesh.i"
 //%import "ncmesh.i"
%import "communication.i"

%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%include "mesh/pncmesh.hpp"

/*
  void Dump(std::ostream &os) const;
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_STDOUT_FILE(ElementSet, Dump)
#endif
