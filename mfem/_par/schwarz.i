%module(package="mfem._par") schwarz
%{
#include "mfem.hpp"      
#include "../../headers/schwarz.hpp"
#include "pyoperator.hpp"    
#include "../common/pycoefficient.hpp"  
#include "numpy/arrayobject.h"    
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array();
%}

%inline %{
#include "../../headers/schwarz.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "pgridfunc.i"
%import "hypre.i"
%import "complex_operator.i"
%import "pmesh.i"

%include "../../headers/schwarz.hpp"

