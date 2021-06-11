/*

   densemat.i

*/
%module (package="mfem._ser") densemat

%feature("autodoc", "1");
%{
#include <fstream>
#include <iostream>
#include "general/zstr.hpp"
#include "linalg/sparsemat.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/io_stream.hpp"
using namespace mfem;  
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%include "exception.i"
%import "mem_manager.i"

%import "array.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "../common/ignore_common_functions.i"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%ignore mfem::DenseMatrix::operator=;
%ignore mfem::DenseTensor::operator=;

%pythonprepend mfem::DenseMatrix::Assign %{
from numpy import ndarray, ascontiguousarray
keep_link = False
if len(args) == 1 and isinstance(args[0], ndarray):
        if args[0].dtype != 'float64':
            raise ValueError('Must be float64 array:' + str(args[0].dtype) + ' was given')
        elif args[0].ndim != 2:
            raise ValueError('Ndim must be two') 
        elif args[0].shape[1] != _densemat.DenseMatrix_Size(self):
            raise ValueError('Length does not match')
        else:
            args = (ascontiguousarray(args[0]),)
%}
%pythonappend mfem::DenseMatrix::Assign %{
    return self
%}
%pythonappend mfem::DenseTensor::Assign %{
    return self
%}
%feature("shadow") mfem::DenseMatrix::operator+= %{
def __iadd__(self, v):
    ret = _densemat.DenseMatrix___iadd__(self, v)
    ret.thisown = self.thisown
    self.thisown = 0                  
    return ret
%}
%feature("shadow") mfem::DenseMatrix::operator-= %{
def __isub__(self, v):
    ret = _densemat.DenseMatrix___isub__(self, v)  
    ret.thisown = self.thisown
    self.thisown = 0            
    return ret
%} 
%feature("shadow") mfem::DenseMatrix::operator*= %{
def __imul__(self, v):
    ret = _densemat.DenseMatrix___imul__(self, v)  
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
  void Assign(PyObject* numpymat) {
    /* note that these error does not raise error in python
       type check is actually done in wrapper layer */
    if (!PyArray_Check(numpymat)){
       PyErr_SetString(PyExc_ValueError, "Input data must be ndarray");
       return;
    }
    int typ = PyArray_TYPE(numpymat);
    if (typ != NPY_DOUBLE){
        PyErr_SetString(PyExc_ValueError, "Input data must be float64");
	return;
    }
    int ndim = PyArray_NDIM(numpymat);
    if (ndim != 2){
      PyErr_SetString(PyExc_ValueError, "Input data NDIM must be two");
      return ;
    }
    npy_intp *shape = PyArray_DIMS(numpymat);    
    int len = self->Width()*self->Height();
    if (shape[1]*shape[0] != len){    
      PyErr_SetString(PyExc_ValueError, "input data length does not match");
      return ;
    }
    PyArrayObject * tmp = 
      PyArray_GETCONTIGUOUS((PyArrayObject *)PyArray_Transpose((PyArrayObject *)numpymat, NULL));
    (* self) = (double *) PyArray_DATA(tmp);
    Py_XDECREF(tmp);
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

/*
  virtual void Print(std::ostream &out = mfem::out, int width_ = 4) const;
  virtual void PrintMatlab(std::ostream &out = mfem::out) const;
  virtual void PrintT(std::ostream &out = mfem::out, int width_ = 4) const;
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(DenseMatrix, Print)
OSTREAM_ADD_DEFAULT_FILE(DenseMatrix, PrintT)
OSTREAM_ADD_DEFAULT_FILE(DenseMatrix, PrintMatlab)				
#endif
