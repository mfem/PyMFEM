%module pmesh
%{
#include <mpi.h>
#include "mesh/pmesh.hpp"
#include "general/communication.hpp"  
#include "numpy/arrayobject.h"
#define MFEM_USE_MPI
  //#include "mpi4py/mpi4py.h"  
%}

%init %{
import_array();
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%import mesh.i
%import pncmesh.i
 //
%import communication.i

%immutable face_nbr_elements;
%immutable face_nbr_vertices;
%immutable gtopo;

#define MFEM_USE_MPI  
%include "mesh/pmesh.hpp"
