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
%import "cpointer.i"
%pointer_class(int, intp);
%import sparsemat.i
%import fespace.i
%import pfespace.i

%ignore DiscreteCurl;
%ignore DiscreteGrad;

%import "hypre_int.i"

%inline %{
  int sizeof_HYPRE_Int(){
    return sizeof(HYPRE_Int);
  }
%}

 /*
    HypreParMatrix(MPI_Comm comm, int nrows, HYPRE_Int glob_nrows,
                  HYPRE_Int glob_ncols, int *I, HYPRE_Int *J,
                  double *data, HYPRE_Int *rows, HYPRE_Int *cols);

    allows to use numpy array to call this
 */
%typemap(in) (int *I, HYPRE_Int *J,
              double *data, HYPRE_Int *rows, HYPRE_Int *cols){

  $1 = (int *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,0)));
  $2 = (HYPRE_Int *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,1)));
  $3 = (double *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,2)));
  $4 = (HYPRE_Int *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,3)));
  $5 = (HYPRE_Int *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,4)));  
}  

%typemap(typecheck ) (int *I, HYPRE_Int *J,
                      double *data, HYPRE_Int *rows,
		      HYPRE_Int *cols){
  /* check if list of 5 numpy array or not */
  if (!PyList_Check($input)) $1 = 0;
  else {
     if (PyList_Size($input) == 5){
       $1 = 1;
       if (!PyArray_Check(PyList_GetItem($input,0))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,1))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,2))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,3))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,4))) $1 = 0;
     } else $1 = 0;       
  }
}

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

