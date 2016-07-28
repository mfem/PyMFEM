%module vector

%{
#include "linalg/vector.hpp"
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "numpy/arrayobject.h"    
%}

// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
 //%import "gridfunc.i"

%typemap(in)  (double *_data, int _size){
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $2 = PyList_Size($input);
  $1 = (double *) malloc(($2)*sizeof(double));
  for (i = 0; i < $2; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (PyInt_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
    } else if (PyFloat_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
    } else {
        free($1);      
        PyErr_SetString(PyExc_ValueError, "List items must be integer/float");
        return NULL;
    }
  }
}
%typemap(typecheck) (double *_data, int _size) {
   $1 = PyList_Check($input) ? 1 : 0;
}


%feature("shadow") mfem::Vector::operator+= %{
def __iadd__(self, v):
    ret = _vector.Vector___iadd__(self, v)
    ret.thisown = self.thisown
    self.thisown = 0                  
    return ret
%}
%feature("shadow") mfem::Vector::operator-= %{
def __isub__(self, v):
    ret = _vector.Vector___isub__(self, v)
    ret.thisown = self.thisown
    self.thisown = 0            
    return ret
%} 
%feature("shadow") mfem::Vector::operator*= %{
def __imul__(self, v):
    ret = _vector.Vector___imul__(self, v)
    ret.thisown = self.thisown
    self.thisown = 0            
    return ret
%} 
%feature("shadow") mfem::Vector::operator/= %{
def __idiv__(self, v):
    ret = _vector.Vector___idiv__(self, v)
    ret.thisown = self.thisown
    self.thisown = 0      
    return ret
%}
%rename(Assign) mfem::Vector::operator=;

// these inlines are to rename add/subtract...
%inline %{
void add_vector(const mfem::Vector &v1, const mfem::Vector &v2, mfem::Vector &v){
   add(v1, v2, v);
}
   /// Do v = v1 + alpha * v2.
void add_vector(const mfem::Vector &v1, double alpha, const mfem::Vector &v2, mfem::Vector &v){
   add(v1, alpha, v2, v);
}
   /// z = a * (x + y)
void add_vector(const double a, const mfem::Vector &x, const mfem::Vector &y, mfem::Vector &z){
   add(a, x, y, z);
}
  /// z = a * x + b * y
void add_vector (const double a, const mfem::Vector &x,
		   const double b, const mfem::Vector &y, mfem::Vector &z){
   add(a, x, b, y, z);
}
   /// Do v = v1 - v2.
void subtract_vector(const mfem::Vector &v1, const mfem::Vector &v2, mfem::Vector &v){
   subtract(v1, v2, v);
}
   /// z = a * (x - y)
void subtract_vector(const double a, const mfem::Vector &x,
		       const mfem::Vector &y, mfem::Vector &z){
   subtract(a, x, y, z);
}
%}

%include "linalg/vector.hpp"

%extend mfem::Vector {
   Vector(const mfem::Vector &v, int offset, int size){
      mfem::Vector *vec;
      vec = new mfem::Vector(v.GetData() +  offset, size);     
      return vec;
  }
  void __setitem__(int i, const double v) {
    (* self)(i) = v;
    }
  const double __getitem__(const int i) const{
    return (* self)(i);
  }
  PyObject* GetDataArray(void) const{
     double * A = self->GetData();    
     int L = self->Size();
     npy_intp dims[] = {L};
     return  PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, A);
  }
};



