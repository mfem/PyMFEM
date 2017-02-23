%module pnonlinearform
%{
#include <mpi.h>
#include "config/config.hpp"
#include "fem/pnonlinearform.hpp"
#include "fem/linearform.hpp"
#include "pyoperator.hpp"           
%}

%include "config/_config.hpp" // include mfem MACRO
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
/*
%init %{
import_array();
%}
*/
%import nonlinearform.i
%import pfespace.i
%import pgridfunc.i

%pointer_class(int, intp);


%include "fem/pnonlinearform.hpp"
