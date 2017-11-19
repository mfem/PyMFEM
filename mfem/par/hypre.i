%module hypre
%{
#include <mpi.h>
#include <Python.h>
#include "fem/gridfunc.hpp"
#include "fem/linearform.hpp"
#include "fem/pfespace.hpp"    
#include "config/config.hpp"        
#include "linalg/hypre.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "iostream_typemap.hpp"    
%}
%include  "config/_config.hpp" // include mfem MACRO

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}
%include "../common/cpointers.i"

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
%typemap(in) (double *_data,  HYPRE_Int *col)(PyArrayObject * tmp_arr1_ = NULL,  PyArrayObject * tmp_arr2_ = NULL){
  //PyArrayObject *tmp_arr1, *tmp_arr2;
  tmp_arr1_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,0));
  tmp_arr2_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,1));
  $1 = (double *) PyArray_DATA(tmp_arr1_);
  $2 = (HYPRE_Int *) PyArray_DATA(tmp_arr2_);
}
%typemap(freearg) (double *_data,  HYPRE_Int *col){
  //Py_XDECREF(tmp_arr1_$argnum); Dont do this.. HypreParVec constructer requires outside object alive
  //Py_XDECREF(tmp_arr2_$argnum);  
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
%typemap(in) (int *I,
	      HYPRE_Int *J,
              double *data,
	      HYPRE_Int *rows,
	      HYPRE_Int *cols)
             (PyArrayObject *tmp_arr1_ = NULL,
	      PyArrayObject *tmp_arr2_ = NULL,
	      PyArrayObject *tmp_arr3_ = NULL,
	      PyArrayObject *tmp_arr4_ = NULL,
	      PyArrayObject *tmp_arr5_ = NULL){
  tmp_arr1_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,0));
  tmp_arr2_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,1));
  tmp_arr3_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,2));
  tmp_arr4_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,3));
  tmp_arr5_ = PyArray_GETCONTIGUOUS((PyArrayObject *)PyList_GetItem($input,4));

  $1 = (int *) PyArray_DATA(tmp_arr1_);
  $2 = (HYPRE_Int *) PyArray_DATA(tmp_arr2_);
  $3 = (double *) PyArray_DATA(tmp_arr3_);
  $4 = (HYPRE_Int *) PyArray_DATA(tmp_arr4_);
  $5 = (HYPRE_Int *) PyArray_DATA(tmp_arr5_);
}
%typemap(freearg) (int *I, HYPRE_Int *J,
		   double *data, HYPRE_Int *rows, HYPRE_Int *cols){
  Py_XDECREF(tmp_arr1_$argnum);
  Py_XDECREF(tmp_arr2_$argnum);  
  Py_XDECREF(tmp_arr3_$argnum);
  Py_XDECREF(tmp_arr4_$argnum);
  Py_XDECREF(tmp_arr5_$argnum);    
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


%newobject mfem::HypreParVector::GlobalVector;
%newobject mfem::HypreParMatrix::Transpose;
%newobject mfem::Add;
%rename(add_hypre) mfem::Add;
%newobject mfem::ParMult;
%newobject mfem::RAP;

%include "linalg/hypre.hpp"

%pythoncode %{
def parvec__repr__(self):
    return "HypreParVector ("+str(self.GlobalSize())+")"
def parvec__del__(self):
    if hasattr(self, "_linked_array"):
        self._linked_arry = None
def parmat__repr__(self):
    shape = (self.GetGlobalNumRows(), self.GetGlobalNumCols())
    lshape = (self.GetNumRows(), self.GetNumCols())  	       
    return "HypreParMatrix "+str(shape)+"["+str(lshape)+"]"
      

HypreParVector.__repr__ = parvec__repr__
HypreParVector.__del__  = parvec__del__      
HypreParMatrix.__repr__ = parmat__repr__
%}
    
%extend mfem::HypreParVector {
PyObject* GetPartitioningArray()
{
  // assumed partitioning mode only
  npy_intp dims[] = {3};
  int typenum =  (sizeof(HYPRE_Int) == 4) ? NPY_INT32 : NPY_INT64;
  HYPRE_Int *part_out;
  
  HYPRE_Int *part = self -> Partitioning();
  PyObject *tmp_arr = PyArray_ZEROS(1, dims, typenum, 0);
  PyObject *arr1 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)tmp_arr);
  Py_XDECREF(tmp_arr);

  part_out = (HYPRE_Int *) PyArray_DATA(arr1);
  part_out[0] = part[0];
  part_out[1] = part[1];
  part_out[2] = part[2];  

  return arr1;
}
}   
%extend mfem::HypreParMatrix {
PyObject* GetRowPartArray()
{
  // assumed partitioning mode only
  npy_intp dims[] = {3};
  int typenum =  (sizeof(HYPRE_Int) == 4) ? NPY_INT32 : NPY_INT64;
  HYPRE_Int *part_out;
  
  HYPRE_Int *part = self -> RowPart();
  PyObject *tmp_arr = PyArray_ZEROS(1, dims, typenum, 0);
  PyObject *arr1 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)tmp_arr);
  Py_XDECREF(tmp_arr);

  part_out = (HYPRE_Int *) PyArray_DATA(arr1);
  part_out[0] = part[0];
  part_out[1] = part[1];
  part_out[2] = part[2];  

  return arr1;
}
PyObject* GetColPartArray()
{
  // assumed partitioning mode only
  npy_intp dims[] = {3};
  int typenum =  (sizeof(HYPRE_Int) == 4) ? NPY_INT32 : NPY_INT64;
  HYPRE_Int *part_out;
  
  HYPRE_Int *part = self -> ColPart();
  PyObject *tmp_arr = PyArray_ZEROS(1, dims, typenum, 0);
  PyObject *arr1 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)tmp_arr);
  Py_XDECREF(tmp_arr);

  part_out = (HYPRE_Int *) PyArray_DATA(arr1);
  part_out[0] = part[0];
  part_out[1] = part[1];
  part_out[2] = part[2];  

  return arr1;
}
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
PyObject* get_local_true_nnz()
{
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
   HYPRE_Int         myid, num_procs, i, j;
   HYPRE_Int         num_nonzeros_offd;
   PyObject *o = NULL;      
   HYPRE_Int tnnz = 0;
   HYPRE_Int nnz;
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
   for (i = 0; i < num_rows; i++)
   {
      for (j = diag_i[i]; j < diag_i[i+1]; j++)
      {
         if ( diag_data )
	 {
           if ((double)diag_data[j] != 0){
    	       tnnz = tnnz + 1;		     
           }
         }		     
		     
      }
      if ( num_nonzeros_offd )
      {
         for (j = offd_i[i]; j < offd_i[i+1]; j++)
         {
            if ( offd_data )
            {
	      if ((double)offd_data[j] != 0){
         	       tnnz = tnnz + 1;		     
               }
     	    }
         }
      }
   }
   o = PyList_New(2);
   PyList_SetItem(o, 0, PyLong_FromLong((long)nnz));
   PyList_SetItem(o, 1, PyLong_FromLong((long)tnnz));

   return o;
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

   PyObject *arr1 = NULL, *arr2 = NULL, *arr3 = NULL, *o = NULL;   
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
   
   HYPRE_Int tnnz;
   tnnz = 0;
   for (i = 0; i < num_rows; i++)
   {
      for (j = diag_i[i]; j < diag_i[i+1]; j++)
      {
         if ( diag_data )
	 {
           if ((double)diag_data[j] != 0){
    	       tnnz = tnnz + 1;		     
           }
         }		     
		     
      }
      if ( num_nonzeros_offd )
      {
         for (j = offd_i[i]; j < offd_i[i+1]; j++)
         {
            if ( offd_data )
            {
	      if ((double)offd_data[j] != 0){
         	       tnnz = tnnz + 1;		     
               }
     	    }
         }
      }
   }
   
   HYPRE_Int* irn;
   HYPRE_Int* jcn;
   double* a;
   //npy_intp dims[] = {nnz};
   npy_intp dims[] = {tnnz};
   int typenum =  (sizeof(HYPRE_Int) == 4) ? NPY_INT32 : NPY_INT64;

   PyObject *tmp_arr1 = PyArray_ZEROS(1, dims, typenum, 0);
   PyObject *tmp_arr2 = PyArray_ZEROS(1, dims, typenum, 0);
   PyObject *tmp_arr3 = PyArray_ZEROS(1, dims, NPY_DOUBLE, 0);
   if (tmp_arr1 == NULL) goto alloc_fail;
   if (tmp_arr2 == NULL) goto alloc_fail;
   if (tmp_arr3 == NULL) goto alloc_fail;

   arr1 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)tmp_arr1);
   arr2 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)tmp_arr2);
   arr3 =  (PyObject *)PyArray_GETCONTIGUOUS((PyArrayObject *)tmp_arr3);
   Py_XDECREF(tmp_arr1);
   Py_XDECREF(tmp_arr2);
   Py_XDECREF(tmp_arr3);
   
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
	   if ((double)diag_data[j] != 0.){
               irn[innz]= I;
               jcn[innz] = J;
               a[innz] = (double)diag_data[j];
               innz = innz + 1;
	   }
         }		     
		     
      }
      if ( num_nonzeros_offd )
      {
         for (j = offd_i[i]; j < offd_i[i+1]; j++)
         {
            J = col_map_offd[offd_j[j]] + base_j;
            if ( offd_data )
            {
     	       if ((double)offd_data[j] != 0.){	      
                   irn[innz]= I;;
                   jcn[innz] = J;
                   a[innz] = (double)offd_data[j];
       	           innz = innz + 1;
	       }
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
     Py_XDECREF(o);     
     return Py_None;
}


}  

