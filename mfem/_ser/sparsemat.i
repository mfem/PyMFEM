%module (package="mfem._ser") sparsemat
%{
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include "linalg/sparsemat.hpp"
#include "numpy/arrayobject.h"
#include "iostream_typemap.hpp"  
#include "pyoperator.hpp"
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%include "exception.i"
%import "array.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "densemat.i"
%import "ostream_typemap.i"
%import "../common/ignore_common_functions.i"
%import "../common/exception.i"

%ignore Walk;
%pythonappend mfem::SparseMatrix::operator*= %{
    val.thisown = 0
    return self
%}
%pythonappend mfem::SparseMatrix::operator+= %{
    val.thisown = 0
    return self
%}
%pythonprepend mfem::SparseMatrix::SparseMatrix %{
import numpy as np  
from scipy.sparse import csr_matrix
if len(args) == 1 and isinstance(args[0], csr_matrix):
   csr = args[0]
   if np.real(csr).dtype != 'float64':
       csr = csr.astype('float64')
   i = np.ascontiguousarray(csr.indptr)
   j = np.ascontiguousarray(csr.indices)
   data = np.ascontiguousarray(csr.data)
   m, n = csr.shape
   this = _sparsemat.new_SparseMatrix([i, j, data, m, n])
   try:
       self.this.append(this)
   except __builtin__.Exception:
       self.this = this
   _sparsemat.SparseMatrix_SetGraphOwner(self, False)
   _sparsemat.SparseMatrix_SetDataOwner(self, False)
   self._i_data = i
   self._j_data = j
   self._d_data = data
  
   return
%}	  

// RAP_P, RAP_R replaces RAP, since RAP has two definition one accept
// pointer and the other accept reference. From Python, two
// can not be distingished..
%inline %{
  mfem::SparseMatrix *RAP_P (const mfem::SparseMatrix &A,
			     const mfem::SparseMatrix &R,
                             mfem::SparseMatrix *ORAP){
    return RAP(A, R, ORAP);
  }    
    
  mfem::SparseMatrix *RAP_R(const mfem::SparseMatrix &Rt,
			    const mfem::SparseMatrix &A,
                            const mfem::SparseMatrix &P){
    return RAP(Rt, A, P);
  }
%}

/*
    support numpy array input to 
    SparseMatrix(int *i, int *j, double *data, int m, int n);
    allows to use numpy array to call this
 */
%typemap(in) (int *i, int *j, double *data, int m, int n)
             (PyArrayObject *tmp_arr1_ = NULL,
	      PyArrayObject *tmp_arr2_ = NULL,
	      PyArrayObject *tmp_arr3_ = NULL,
	      int tmp_4_, int tmp_5_ ){
  tmp_arr1_ = (PyArrayObject *)PyList_GetItem($input,0);
  tmp_arr2_ = (PyArrayObject *)PyList_GetItem($input,1);
  tmp_arr3_ = (PyArrayObject *)PyList_GetItem($input,2);
  tmp_4_ = PyInt_AsLong(PyList_GetItem($input,3));
  tmp_5_ = PyInt_AsLong(PyList_GetItem($input,4));

  $1 = (int *) PyArray_DATA(tmp_arr1_);
  $2 = (int *) PyArray_DATA(tmp_arr2_);
  $3 = (double *) PyArray_DATA(tmp_arr3_);
  $4 = (int) tmp_4_;
  $5 = (int) tmp_5_;
  //PyArray_CLEARFLAGS(tmp_arr1_, NPY_ARRAY_OWNDATA);
  //PyArray_CLEARFLAGS(tmp_arr2_, NPY_ARRAY_OWNDATA);
  //
  PyArray_CLEARFLAGS(tmp_arr3_, NPY_ARRAY_OWNDATA);    
}
%typemap(freearg) (int *i, int *j, double *data, int m, int n){
  //Py_XDECREF(tmp_arr1_$argnum); Dont do this.. We set OwnsGraph and OwnsData to Fase in Python
  //Py_XDECREF(tmp_arr2_$argnum);  
  //Py_XDECREF(tmp_arr3_$argnum);
}

%typemap(typecheck ) (int *i, int *j, double *data, int m, int n){
  /* check if list of 5 numpy array or not */
  if (!PyList_Check($input)) $1 = 0;
  else {
     if (PyList_Size($input) == 5){
       $1 = 1;
       if (!PyArray_Check(PyList_GetItem($input,0))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,1))) $1 = 0;
       if (!PyArray_Check(PyList_GetItem($input,2))) $1 = 0;
       if (!PyInt_Check(PyList_GetItem($input,3))) $1 = 0;
       if (!PyInt_Check(PyList_GetItem($input,4))) $1 = 0;
     } else $1 = 0;       
  }
}


%include "linalg/sparsemat.hpp"

%extend mfem::SparseMatrix {
PyObject* GetIArray(void) const{
  int * I = self->GetI();
  int L = self->Size();

  npy_intp dims[] = { L+1 };
  return  PyArray_SimpleNewFromData(1, dims, NPY_INT, I);
  }
PyObject* GetJArray(void) const{
  int * I = self->GetI();
  int * J = self->GetJ();
  int L = self->Size();
  npy_intp dims[] = { I[L]};
  return  PyArray_SimpleNewFromData(1, dims, NPY_INT, J);
  }
PyObject* GetDataArray(void) const{
  int * I = self->GetI();
  double * A = self->GetData();    
  int L = self->Size();
  npy_intp dims[] = {I[L]};
  return  PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, A);
  }
void Print(const char *file){
  std::ofstream ofile(file);
  if (!ofile)
  {
     std::cerr << "\nCan not produce output file: " << file << '\n' << std::endl;
     return;
  }
  self -> Print(ofile);
  ofile.close();
  }
};
