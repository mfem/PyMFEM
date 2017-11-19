%module array
%rename(Equal) mfem::Array <class T>::operator=;
%{
#include <iostream>
#include <stdio.h>
#include "iostream_typemap.hpp"    
#include "general/array.hpp"  
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%import "ostream_typemap.i"

%typemap(in) (int *_data, int asize) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $2 = PyList_Size($input);
  $1 = (int *) malloc(($2)*sizeof(int));
  for (i = 0; i < $2; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyInt_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be integer");
        return NULL;
    }
    $1[i] = (int)PyInt_AsLong(s);
  }
}
%typemap(typecheck) (int *_data, int asize) {
   $1 = PyList_Check($input) ? 1 : 0;
}

%typemap(newfree) (int *_data,  int asize) {
   if ($1) free($1);
}

%pythonappend mfem::Array::Array %{
  if len(args) == 1 and isinstance(args[0], list):
      self.MakeDataOwner()
%}

//%import "intrules.i"
//%newobject intArray
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

};
namespace mfem{
%feature("shadow")Array::FakeToList %{
def ToList(self):
    return [self[i] for i in range(self.Size())]
%}
}  
namespace mfem{
%template(intArray) Array<int>;
%template(doubleArray) Array<double>;
%template(doubleSwap) Swap<double>;  
%template(intSwap) Swap<int>;  
}




