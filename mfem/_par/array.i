%module(package="mfem._par") array
%feature("autodoc", "1");

//%rename(Equal) mfem::Array <class T>::operator=;
%{
#include <fstream>  
#include <iostream>
#include <stdio.h>
#include <vector>  
#include "../common/io_stream.hpp"
#include "mfem.hpp"
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
%import "../common/memorytype_typemap.i"
ENUM_TO_MEMORYTYPE(mfem::MemoryType mt)

 /*  boolArray, intArray, doubleArray */
%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT(int, PyLong_AsLong)
ARRAY_LISTTUPLE_INPUT(double, PyFloat_AsDouble)
ARRAY_LISTTUPLE_INPUT(bool, PyObject_IsTrue)

%import "../common/data_size_typemap.i"
XXXPTR_SIZE_IN(int *data_, int asize, int)
XXXPTR_SIZE_IN(double *data_, int asize, double)
XXXPTR_SIZE_IN(bool *data_, int asize, bool)

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
  mfem::Array (void *List_or_Tuple, T *_unused){
    /*
    This method is wrapped to recived tuple or list to create
    Array object
    */
    mfem::Array <T>  *arr;
    arr = new mfem::Array<T>(*(int*)List_or_Tuple);    
    return arr;
  }
  void __setitem__(int i, const T v) {
    (* self)[i] = v;
    }
  void Assign(const T &a){
     *self = a;
  }   
  void FakeToList(void){}
  void __iter__(void){}
};
namespace mfem{
%pythonprepend Array::__setitem__ %{
    i = int(i)
%}
  //%pythonprepend Array::__getitem__ %{
  //    i = int(i)
  //%}  
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
%template(doubleSwap) Swap<double>;  
%template(intSwap) Swap<int>;  
}
%import "../common/array_instantiation_macro.i"
INSTANTIATE_ARRAY_INT
INSTANTIATE_ARRAY_DOUBLE
IGNORE_ARRAY_METHODS(bool)
INSTANTIATE_ARRAY_BOOL


