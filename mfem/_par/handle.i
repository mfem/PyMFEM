%module(package="mfem._par") handle
%feature("autodoc", "1");
%{
#include <fstream>  
#include <iostream>
#include "../common/io_stream.hpp"        
#include "config/config.hpp"  
#include "linalg/hypre.hpp"
#include "linalg/handle.hpp"  
#include "fem/gridfunc.hpp"  
#include "fem/linearform.hpp"
#include "fem/pfespace.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
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
%import "../common/exception.i"

%import "operators.i"

#ifdef MFEM_USE_MPI
%import "hypre.i"
#endif
#ifdef MFEM_USE_PETSC
%include "petsc.i"
#endif

%import "mem_manager.i"

%import "../common/handle_template.i"

// instatitate template methods (step 1: Rename Macro )

 //CONSTRUCTOR_RENAME(SparseMatrix)
AS_RENAME(SparseMatrix)
IS_RENAME(SparseMatrix)
GET_RENAME(SparseMatrix)
RESET_RENAME(SparseMatrix)
CONVERT_FROM_RENAME(SparseMatrix)

 //CONSTRUCTOR_RENAME(HypreParMatrix)  
AS_RENAME(HypreParMatrix)
IS_RENAME(HypreParMatrix)
GET_RENAME(HypreParMatrix)
RESET_RENAME(HypreParMatrix)
CONVERT_FROM_RENAME(HypreParMatrix)

//CONSTRUCTOR_RENAME(ComplexHypreParMatrix)  
AS_RENAME(ComplexHypreParMatrix)
IS_RENAME(ComplexHypreParMatrix)
GET_RENAME(ComplexHypreParMatrix)
RESET_RENAME(ComplexHypreParMatrix)
CONVERT_FROM_RENAME(ComplexHypreParMatrix)

#ifdef MFEM_USE_PETSC
 //CONSTRUCTOR_RENAME(PetscParMatrix)  
AS_RENAME(PetscParMatrix)
IS_RENAME(PetscParMatrix)
GET_RENAME(PetscParMatrix)
RESET_RENAME(PetscParMatrix)
CONVERT_FROM_RENAME(PetscParMatrix)  
#endif

%include "linalg/handle.hpp"

%pythoncode %{
OperatorPtr=OperatorHandle  
%}

// instatitate template methods (step 2: Instantiation)
//CONSTRUCTOR_WRAP(SparseMatrix)
AS_WRAP(SparseMatrix)
IS_WRAP(SparseMatrix)
GET_WRAP(SparseMatrix)
RESET_WRAP(SparseMatrix)          
CONVERT_FROM_WRAP(SparseMatrix)

//CONSTRUCTOR_WRAP(HypreParMatrix)  
AS_WRAP(HypreParMatrix)
IS_WRAP(HypreParMatrix)
GET_WRAP(HypreParMatrix)
RESET_WRAP(HypreParMatrix)          
CONVERT_FROM_WRAP(HypreParMatrix)

//CONSTRUCTOR_WRAP(ComplexHypreParMatrix)  
AS_WRAP(ComplexHypreParMatrix)
IS_WRAP(ComplexHypreParMatrix)
GET_WRAP(ComplexHypreParMatrix)
RESET_WRAP(ComplexHypreParMatrix)          
CONVERT_FROM_WRAP(ComplexHypreParMatrix)
  
#ifdef MFEM_USE_PETSC
//CONSTRUCTOR_WRAP(PetscParMatrix)    
AS_WRAP(PetscParMatrix)
IS_WRAP(PetscParMatrix)
GET_WRAP(PetscParMatrix)
RESET_WRAP(PetscParMatrix)          
CONVERT_FROM_WRAP(PetscParMatrix)
#endif
  
