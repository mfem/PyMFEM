%module pncmesh
%{
#include <mpi.h>
#include "config/config.hpp"
#include "mesh/mesh_headers.hpp"
#include "mpi4py/mpi4py.h"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"  
%}

%include "config/_config.hpp" // include mfem MACRO
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}


%init %{
import_array();
%}

%import mesh.i
%import ncmesh.i
%import communication.i

%pointer_class(int, intp);

%include "mesh/pncmesh.hpp"
