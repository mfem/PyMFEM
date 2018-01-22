%module pmesh
%{
#include <mpi.h>
#include "iostream_typemap.hpp"     
#include "config/config.hpp"
#include "mesh/pmesh.hpp"
#include "fem/linearform.hpp"
#include "general/communication.hpp"  
#include "numpy/arrayobject.h"
%}

%include  "config/_config.hpp" // include mfem MACRO
%init %{
import_array();
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%include "../common/cpointers.i"
%import mesh.i
%import pncmesh.i
 //
%import communication.i
%import "ostream_typemap.i"

%immutable face_nbr_elements;
%immutable face_nbr_vertices;
%immutable gtopo;

%pointer_class(int, intp);

%include "mesh/pmesh.hpp"
