%module(package="mfem._par", directors="1") tmop
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%include "exception.i"
%import "../common/exception_director.i"
%include "../common/typemap_macros.i"
LIST_TO_MFEMOBJ_BOOLARRAY_IN(const mfem::Array<bool>& )

%import intrules.i
%import gridfunc.i

%import "../common/array_instantiation_macro.i"
IGNORE_ARRAY_METHODS(mfem::TMOP_Integrator *)
INSTANTIATE_ARRAY0(TMOP_Integrator *, TMOP_Integrator, 1) 

/* destructor handles freeing evaluator */
%pythonappend mfem::DiscreteAdaptTC::SetAdaptivityEvaluator %{
    ae.thisown = 0
%}

%include "fem/tmop.hpp"

