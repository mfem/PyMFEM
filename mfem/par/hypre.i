%module hypre
%{
#include <mpi.h>
#include <Python.h>
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
%import vector.i
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
   support numpy array input to 
   HypreParVector(MPI_Comm comm, HYPRE_Int glob_size, double *_data,
                  HYPRE_Int *col);

*/
%typemap(in) (double *_data,  HYPRE_Int *col){
  $1 = (double *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,0)));
  $2 = (HYPRE_Int *) PyArray_DATA(PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,1)));
}
%typemap(typecheck )(double *_data,  HYPRE_Int *col){
  /* check if list of 2 numpy array or not */
  if (!PyList_Check($input)) $1 = 0;
  else {
     if (PyList_Size($input) == 2){
       $1 = 1;
       if (!PyArray_Check(PyList_GetItem($input,0))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,1))) $1 = 0;
     } else $1 = 0;       
  }
}
 /*
    support numpy array input to 
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
%pythonappend mfem::HypreParVector::HypreParVector %{
  if isinstance(args[-1], list):
     # in this case, ParVector does not own the object
     # in order to prevent python from freeing the input
     # array, object is kept in ParVector
     # args[-1][0]  _data
     # args[-1][0]  col
     self._linked_array = args[-1][0]
%}
%pythonappend mfem::HypreParMatrix::operator*= %{
#    val.thisown = 0
    return self
%}
%include "linalg/hypre.hpp"

%extend mfem::HypreParVector {
PyObject* GetPartitioningArray()
{
  // assumed partitioning mode only
  npy_intp dims[] = {3};
  int typenum =  (sizeof(HYPRE_Int) == 4) ? NPY_INT32 : NPY_INT64;
  HYPRE_Int *part_out;
  
  HYPRE_Int *part = self -> Partitioning();
  PyObject *arr1 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)PyArray_ZEROS(1, dims, typenum, 0));

  part_out = (HYPRE_Int *) PyArray_DATA(arr1);
  part_out[0] = part[0];
  part_out[1] = part[1];
  part_out[2] = part[2];  

  return arr1;
}
}   
%extend mfem::HypreParMatrix {
HYPRE_Int get_local_nnz()//mfem::HypreParMatrix *pmatrix)
{
  //hypre_ParCSRMatrix *matrix =  static_cast<hypre_ParCSRMatrix *>(*pmatrix);
   hypre_ParCSRMatrix *matrix =  static_cast<hypre_ParCSRMatrix *>(*self);
   MPI_Comm          comm;
   hypre_CSRMatrix  *diag;
   hypre_CSRMatrix  *offd;
   if (!matrix)
   {
      /*hypre_error_in_arg(1);*/
     return 0;
   }
   comm = hypre_ParCSRMatrixComm(matrix);
   diag            = hypre_ParCSRMatrixDiag(matrix);
   offd            = hypre_ParCSRMatrixOffd(matrix);
   return hypre_CSRMatrixNumNonzeros(diag) + hypre_CSRMatrixNumNonzeros(offd);
}
 
/* MakeMatrixCoordinateFormat */
/*   generated from HYPRE_TO_MUMPS, but returns 0 based array */
PyObject* GetCooDataArray(const HYPRE_Int           base_i = 0,
			  const HYPRE_Int           base_j = 0)
{
  //hypre_ParCSRMatrix *matrix =  static_cast<hypre_ParCSRMatrix *>(*pmatrix);
   hypre_ParCSRMatrix *matrix =  static_cast<hypre_ParCSRMatrix *>(*self);


   MPI_Comm          comm;
   HYPRE_Int         first_row_index;
   HYPRE_Int         first_col_diag;
   hypre_CSRMatrix  *diag;
   hypre_CSRMatrix  *offd;
   HYPRE_Int        *col_map_offd;
   HYPRE_Int         num_rows;
   HYPRE_Int        *row_starts;
   HYPRE_Int        *col_starts;
   HYPRE_Complex    *diag_data;
   HYPRE_Int        *diag_i;
   HYPRE_Int        *diag_j;
   HYPRE_Complex    *offd_data;
   HYPRE_Int        *offd_i;
   HYPRE_Int        *offd_j;
   HYPRE_Int         myid, num_procs, i, j, I, J;
   HYPRE_Int         num_nonzeros_offd;
   HYPRE_Int         ilower, iupper, jlower, jupper;
   
   HYPRE_Int innz = 0;
   HYPRE_Int nnz;

   PyObject *arr1, *arr2, *arr3, *o;   
   if (!matrix)
   {
      /*hypre_error_in_arg(1);*/
     return Py_None;    
   }

   comm = hypre_ParCSRMatrixComm(matrix);
   diag = hypre_ParCSRMatrixDiag(matrix);
   offd = hypre_ParCSRMatrixOffd(matrix);
   nnz =  hypre_CSRMatrixNumNonzeros(diag) + hypre_CSRMatrixNumNonzeros(offd);

   first_row_index = hypre_ParCSRMatrixFirstRowIndex(matrix);
   first_col_diag  = hypre_ParCSRMatrixFirstColDiag(matrix);
   diag            = hypre_ParCSRMatrixDiag(matrix);
   offd            = hypre_ParCSRMatrixOffd(matrix);
   col_map_offd    = hypre_ParCSRMatrixColMapOffd(matrix);
   num_rows        = hypre_ParCSRMatrixNumRows(matrix);
   row_starts      = hypre_ParCSRMatrixRowStarts(matrix);
   col_starts      = hypre_ParCSRMatrixColStarts(matrix);
   hypre_MPI_Comm_rank(comm, &myid);
   hypre_MPI_Comm_size(comm, &num_procs);
   num_nonzeros_offd = hypre_CSRMatrixNumNonzeros(offd);

   diag_data = hypre_CSRMatrixData(diag);
   diag_i    = hypre_CSRMatrixI(diag);
   diag_j    = hypre_CSRMatrixJ(diag);
   offd_i    = hypre_CSRMatrixI(offd);
   if (num_nonzeros_offd)
   {
      offd_data = hypre_CSRMatrixData(offd);
      offd_j    = hypre_CSRMatrixJ(offd);
   }

   //#ifdef HYPRE_NO_GLOBAL_PARTITION
   if (HYPRE_AssumedPartitionCheck()){
     //std::cout << "no_global_partition\n";
      ilower = row_starts[0]+base_i;
      iupper = row_starts[1]+base_i - 1;
      jlower = col_starts[0]+base_j;
      jupper = col_starts[1]+base_j - 1;
   } else {
      ilower = row_starts[myid]  +base_i;
      iupper = row_starts[myid+1]+base_i - 1;
      jlower = col_starts[myid]  +base_j;
      jupper = col_starts[myid+1]+base_j - 1;
   }
   
   HYPRE_Int* irn;// = (HYPRE_Int *)malloc(sizeof(HYPRE_Int)*nnz);
   HYPRE_Int* jcn;// = (HYPRE_Int *)malloc(sizeof(HYPRE_Int)*nnz);
   double* a;// = (double *)malloc(sizeof(double)*nnz);
   npy_intp dims[] = {nnz};
   int typenum =  (sizeof(HYPRE_Int) == 4) ? NPY_INT32 : NPY_INT64;

   //std::cout << "nnz " << std::to_string(nnz) << "\n";
   arr1 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)PyArray_ZEROS(1, dims, typenum, 0));
   arr2 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)PyArray_ZEROS(1, dims, typenum, 0));
   arr3 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)PyArray_ZEROS(1, dims, NPY_DOUBLE, 0));
   
   if (arr1 == NULL) goto alloc_fail;
   if (arr2 == NULL) goto alloc_fail;
   if (arr3 == NULL) goto alloc_fail;
   irn = (HYPRE_Int *) PyArray_DATA(arr1);
   jcn = (HYPRE_Int *) PyArray_DATA(arr2);
   a = (double *) PyArray_DATA(arr3);
   //*irn_p = irn;
   //*jcn_p = jcn;
   //*a_p   = a;

   for (i = 0; i < num_rows; i++)
   {
      I = first_row_index + i + base_i;
      for (j = diag_i[i]; j < diag_i[i+1]; j++)
      {
         J = first_col_diag + diag_j[j] + base_j;
         if ( diag_data )
	 {
           irn[innz]= I;
           jcn[innz] = J;
           a[innz] = (double)diag_data[j];
	   innz = innz + 1;		     
         }		     
		     
      }
      if ( num_nonzeros_offd )
      {
         for (j = offd_i[i]; j < offd_i[i+1]; j++)
         {
            J = col_map_offd[offd_j[j]] + base_j;
            if ( offd_data )
            {
               irn[innz]= I;;
               jcn[innz] = J;
               a[innz] = (double)offd_data[j];
    	       innz = innz + 1;		     	       
     	    }
         }
      }
   }

  o = PyList_New(8);
  PyList_SetItem(o, 0, PyLong_FromLong((long)num_rows));
  PyList_SetItem(o, 1, PyLong_FromLong((long)ilower));
  PyList_SetItem(o, 2, PyLong_FromLong((long)iupper));
  PyList_SetItem(o, 3, PyLong_FromLong((long)jlower));
  PyList_SetItem(o, 4, PyLong_FromLong((long)jupper));  
  PyList_SetItem(o, 5, arr1);
  PyList_SetItem(o, 6, arr2);
  PyList_SetItem(o, 7, arr3);
  
  return o;
  alloc_fail:
     Py_XDECREF(arr1);
     Py_XDECREF(arr2);
     Py_XDECREF(arr3);
     return Py_None;
}


}  

