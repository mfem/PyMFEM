%module(package="mfem._ser") array

%feature("autodoc", "1");

%rename(Equal) mfem::Array <class T>::operator=;
%{
#include <fstream>  
#include <iostream>
#include <stdio.h>
#include "../common/io_stream.hpp"
#include "general/zstr.hpp"
#include "general/array.hpp"  
#include "numpy/arrayobject.h"    
%}

%begin %{
#define PY_SSIZE_T_CLEAN
%}
%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%import "mem_manager.i"

%import "../common/data_size_typemap.i"
INTPTR_SIZE_IN(int *data_, int asize)
DOUBLEPTR_SIZE_IN(double *data_, int asize)

// intArray constructor
 /*
%typemap(in) (int *data_, int asize) (int * temp_ptr, bool ptr_given=false, PyObject *s1, PyObject *s2){
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $2 = PyList_Size($input);
  if ($2 == 2){
     s1 = PyList_GetItem($input,0);
     s2 = PyList_GetItem($input,1);
     if (SWIG_ConvertPtr(s1, (void **) &temp_ptr, $descriptor(const int *), 0 |0) == -1) {
        ptr_given=false;
     } else {
        ptr_given=true;
     }
  }
  if (! ptr_given){
       $1 = (int *) malloc(($2)*sizeof(int));    
      for (i = 0; i < $2; i++) {
         PyObject *s = PyList_GetItem($input,i);
         if (PyInt_Check(s)) {
             $1[i] = (int)PyInt_AsLong(s);
         } else if ((PyArray_PyIntAsInt(s) != -1) || !PyErr_Occurred()) {
             $1[i] = PyArray_PyIntAsInt(s);
	 } else {
             free($1);
             PyErr_SetString(PyExc_ValueError, "List items must be integer");
             return NULL;
	 }
      }
  } else {
    $1 = temp_ptr;
    $2 = PyLong_AsLong(s2);    
  }
}
%typemap(typecheck) (int *data_, int asize) {
   $1 = PyList_Check($input) ? 1 : 0;
}

%typemap(newfree) (int *data_,  int asize) {
   if ($1) free($1);
}

%typemap(in) (double *data_, int asize) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $2 = PyList_Size($input);
  $1 = (double *) malloc(($2)*sizeof(int));
  for (i = 0; i < $2; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (PyInt_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
    } else if (PyFloat_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
    } else {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be integer");
        return NULL;
    }
  }
}
%typemap(typecheck) (double *data_, int asize) {
   $1 = PyList_Check($input) ? 1 : 0;
}

%typemap(newfree) (double *data_,  int asize) {
   if ($1) free($1);
}
*/

%pythonappend mfem::Array::Array %{
  if len(args) == 1 and isinstance(args[0], list):
      if (len(args[0]) == 2 and hasattr(args[0][0], 'disown') and
	  not hasattr(args[0][1], 'disown')):
          ## first element is SwigObject, like <Swig Object of type 'int *'>
          ## We do not own data in this case.
          pass
      else:
          self.MakeDataOwner()
%}

//%import "intrules.i"
//%newobject intArray
%ignore mfem::Array::operator[];
%ignore mfem::Array2D::operator[];
%ignore mfem::BlockArray::operator[];
%ignore mfem::Array::operator=;
%ignore mfem::Array2D::operator=;
%ignore mfem::Array::operator T *;
%ignore mfem::Array::operator const T *;

%include "general/array.hpp"

%extend mfem::Array{
  void __setitem__(int i, const T v) {
    (* self)[i] = v;
    }
  const T & __getitem__(const int i) const{
    return (* self)[i];
  }
  void Assign(const T &a){
     *self = a;
  }
  void FakeToList(void){}
  /* since Array is template class I can not fill PyList object
     in C++ side */
  void __iter__(void){}
  /* this will be shadowed by Python */
};
namespace mfem{
%pythonprepend Array::__setitem__ %{
    i = int(i)
%}
%pythonprepend Array::__getitem__ %{
    i = int(i)
%}
%feature("shadow")Array::FakeToList %{
def ToList(self):
    return [self[i] for i in range(self.Size())]
%}
%feature("shadow")Array::__iter__ %{
def __iter__(self):
    class iter_array:
        def __init__(self, obj):
            self.obj = obj
            self.idx = 0
            self.size = obj.Size()
        def __iter__(self):
            self.idx = 0
        def __next__(self):
            if self.idx < self.size:
                res = self.obj[self.idx]
                self.idx += 1
                return res
            else:
                raise StopIteration
    return iter_array(self)
%}
}

/*
void Print(std::ostream &out = mfem::out, int width = 4) const;
void Save(std::ostream &out, int fmt = 0) const;
void Save(std::ostream &out, int fmt = 0) const
void Print(std::ostream &out = mfem::out, int width = 4);
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(Array, Print)
OSTREAM_ADD_DEFAULT_FILE(Array2D, Print)
OSTREAM_ADD_DEFAULT_STDOUT_FILE(Array, Save)
OSTREAM_ADD_DEFAULT_STDOUT_FILE(Array2D, Save)
#endif

namespace mfem{
%template(intArray) Array<int>;
%template(doubleArray) Array<double>;
%template(doubleSwap) Swap<double>;  
%template(intSwap) Swap<int>;  
}


