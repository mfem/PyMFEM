%module(package="mfem._par") pnonlinearform
%{
#include <mpi.h>
#include "config/config.hpp"
#include "fem/pnonlinearform.hpp"
#include "fem/linearform.hpp"
#include "numpy/arrayobject.h"  
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
   import_array();
%}

%include "exception.i"
%import "vector.i"
%import "nonlinearform.i"
%import "blockoperator.i"
%import "pfespace.i"
%import "pgridfunc.i"
%import "../common/exception.i"
%include "../common/typemap_macros.i"

%pointer_class(int, intp);

LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::ParFiniteElementSpace *> &pf,
    		        mfem::ParFiniteElementSpace *)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Array<int> *> &bdr_attr_is_ess,
 		        mfem::Array<int> *)
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::Vector *> &rhs, mfem::Vector *)

%include "fem/pnonlinearform.hpp"
