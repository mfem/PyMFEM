%module densemat

%{
#include "linalg/sparsemat.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "iostream_typemap.hpp"
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
//%import "general/array.hpp"
%import "array.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "ostream_typemap.i"
%import "../common/ignore_common_functions.i"

%ignore mfem::DenseMatrix::operator=;
%pythonappend mfem::DenseMatrix::Assign %{
    return self
%}
%pythonappend mfem::DenseTensor::Assign %{
    return self
%}

%feature("shadow") mfem::DenseMatrix::operator+= %{
def __iadd__(self, v):
    ret = _densmat.DenseMatrix___iadd__(self, v)
    ret.thisown = self.thisown
    self.thisown = 0                  
    return ret
%}
%feature("shadow") mfem::DenseMatrix::operator-= %{
def __isub__(self, v):
    ret = _densmat.DenseMatrix___isub__(self, v)  
    ret.thisown = self.thisown
    self.thisown = 0            
    return ret
%} 
%feature("shadow") mfem::DenseMatrix::operator*= %{
def __imul__(self, v):
    ret = _densmat.DenseMatrix___imul__(self, v)  
    ret.thisown = self.thisown
    self.thisown = 0            
    return ret
%} 
%feature("shadow") mfem::DenseMatrix::__setitem__%{
def __setitem__(self, *args):
    i, j, v = args[0][0], args[0][1], args[1]
    return _densemat.DenseMatrix___setitem__(self, i, j, v)
%}
%feature("shadow") mfem::DenseMatrix::__getitem__%{
def __getitem__(self, *args):
    i, j = args[0][0], args[0][1]
    return _densemat.DenseMatrix___getitem__(self, i, j)
%}
%feature("shadow") mfem::DenseTensor::__setitem__%{
def __setitem__(self, *args):
    i, j, k, v = args[0][0], args[0][1], args[0][2], args[1]
    return _densemat.DenseTensor___setitem__(self, i, j, k, v)
%}
%feature("shadow") mfem::DenseTensor::__getitem__%{
def __getitem__(self, *args):
  try:
     check = len(args[0]) == 3
  except:
     check = False
  if check:
     i, j, k = args[0][0], args[0][1], args[0][2]
     return _densemat.DenseTensor___getitem__(self, i, j, k)
  try:
     check = int(args[0])
  except:
     check = -1
  if check >= 0:     
     return _densemat.DenseTensor___getitem__(self, check)
%} 

%include "linalg/densemat.hpp"

%extend mfem::DenseMatrix {
  void Assign(const double v) {
    (* self) = v;
  }
  void Assign(const mfem::DenseMatrix &m) {
    (* self) = m;
  }
  const double __getitem__(const int i, const int j) const{
    return (* self)(i, j);
  }
  void __setitem__(int i, int j,  const double v) {
    (* self)(i, j) = v;
  }
  PyObject* GetDataArray(void) const{
     double * A = self->Data();    
     npy_intp dims[] = {self->Width(), self->Height()};
     return  PyArray_Transpose((PyArrayObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, A), NULL);
  }
};
%extend mfem::DenseTensor {
  void Assign(const double c) {
    (* self) = c;
  }
  const double __getitem__(const int i, const int j, const int k) const{
    return (* self)(i, j, k);
  }
  const mfem::DenseMatrix& __getitem__(const int k) const{
    return (* self)(k);
  }
  void __setitem__(int i, int j, int k, const double v) {
    (* self)(i, j, k) = v;
  }
  PyObject* GetDataArray(void){
     // DoDo this method can not be const since DenseTensor::Data is not const
     double * A = self->Data();    
     npy_intp dims[] = {self->SizeK(), self->SizeJ(), self->SizeI()};
     PyObject * obj = PyArray_SimpleNewFromData(3, dims, NPY_DOUBLE, A);
     //obj = PyArray_SwapAxes((PyArrayObject *)obj, 0, 2);
     obj = PyArray_SwapAxes((PyArrayObject *)obj, 1, 2);
     return obj;
  }
};

