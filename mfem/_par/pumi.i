%module(package="mfem._par", directors="0")  pumi

%feature("autodoc", "1");

%{
#include "mesh/mesh_headers.hpp"
#include "mesh/pumi.hpp"
#include "fem/fem.hpp"
#include "general/array.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "../common/io_stream.hpp"         
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"

%}

%include "../common/mfem_config.i"

#ifndef MFEM_USE_MPI
   #define MFEM_USE_MPI  YES
#endif

#ifndef MFEM_USE_PUMI
   #define MFEM_USE_PUMI YES
#endif

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array();
%}

//%import "cpointers.i"
%include "exception.i"
%import "../common/numpy_int_typemap.i"
%import "pgridfunc.i"
%import "mesh.i"
%import "pmesh.i"

%include "../common/exception.i"

%inline %{
  mfem::ParPumiMesh *ParMesh2ParPumiMesh(mfem::ParMesh *pmesh) {
    return dynamic_cast<mfem::ParPumiMesh*>(pmesh);
  }
%}

#undef MFEM_USE_PUMI
#define MFEM_USE_PUMI YES

%include "mesh/pumi.hpp"



