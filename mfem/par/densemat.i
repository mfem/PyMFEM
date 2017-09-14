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
%import "general/array.hpp"
%import "array.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "ostream_typemap.i"

%rename(Assign) mfem::DenseMatrix::operator=;
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

%rename(add_dense) mfem::Add;
%include "linalg/densemat.hpp"

%extend mfem::DenseMatrix {
  const double __getitem__(const int i, const int j) const{
    return (* self)(i, j);
  }
  void __setitem__(int i, int j,  const double v) {
    (* self)(i, j) = v;
  }
  PyObject* GetDataArray(void) const{
     double * A = self->Data();    
     npy_intp dims[] = {self->Width(), self->Height()};
     return  PyArray_CopyAndTranspose(PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, A));
  }
};

