%module pncmesh
%{
#include <mpi.h>
#include "mesh/mesh_headers.hpp"
#include "mpi4py/mpi4py.h"
#include "numpy/arrayobject.h"  
#define MFEM_USE_MPI  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#define MFEM_USE_MPI
%init %{
import_array();
%}


%init %{
import_array();
%}

%import mesh.i
%import ncmesh.i
%import communication.i
#define MFEM_USE_MPI  
%include "mesh/pncmesh.hpp"
