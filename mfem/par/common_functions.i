%module common_functions

%{
#include "iostream_typemap.hpp"
#include "linalg/blockmatrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/densemat.hpp"
#include "linalg/hypre.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"     
  %}
%include  "config/_config.hpp" // include mfem MACRO (MFEM_USE_MPI)

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}
%import "array.i"

%ignore "";
%rename("%s") mfem;
%rename("%s") mfem::Add;
%rename("%s") mfem::Transpose;
%rename("%s") mfem::Mult;
%rename("%s") mfem::RAP;
%rename("%s") mfem::InnerProduct;
%newobject mfem::Add;
%newobject mfem::RAP;
%include "general/table.hpp"
%include "linalg/vector.hpp"
%include "linalg/blockmatrix.hpp"
%include "linalg/densemat.hpp"
%include "linalg/sparsemat.hpp"
%include "linalg/hypre.hpp"

