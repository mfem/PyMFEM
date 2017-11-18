%module plinearform
%{
#include "config/config.hpp"
#include "fem/gridfunc.hpp"
#include "fem/pgridfunc.hpp"  
#include "fem/plinearform.hpp"
#include "numpy/arrayobject.h"
%}
%include  "config/_config.hpp" // include mfem MACRO

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}

%import linearform.i
%import pfespace.i
%import pgridfunc.i
%import hypre.i

%newobject mfem::ParLinearForm::ParallelAssemble;
%pointer_class(int, intp);
%include "fem/plinearform.hpp"

