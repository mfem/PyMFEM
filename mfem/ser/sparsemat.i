%module sparsemat

%{
#include "linalg/sparsemat.hpp"
#include "numpy/arrayobject.h"  
  %}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "densemat.i"
%ignore Walk;

%pythonappend mfem::SparseMatrix::operator*= %{
    val.thisown = self.thisown
%}
%pythonappend mfem::SparseMatrix::operator+= %{
    val.thisown = self.thisown
%}

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
};
