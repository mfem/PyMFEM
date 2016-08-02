%module hypre
%{
#include <mpi.h>
#include "linalg/hypre.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"           
#define MFEM_USE_MPI

%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}

%import sparsemat.i
%import fespace.i
%import pfespace.i

%ignore DiscreteCurl;
%ignore DiscreteGrad;

%import "hypre_int.i"
/*
%inline %{
#ifdef HYPRE_BIGINT
typedef long long int HYPRE_Int;
#define HYPRE_MPI_INT MPI_LONG_LONG_INT
#else 
typedef int HYPRE_Int;
#define HYPRE_MPI_INT MPI_INT
#endif
%}
*/
#define MFEM_USE_MPI  
%pythonappend mfem::HypreParMatrix::operator*= %{
#    val.thisown = 0
    return self
%}
%include "linalg/hypre.hpp"

