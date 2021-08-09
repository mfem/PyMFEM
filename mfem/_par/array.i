%module(package="mfem._par") array

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
  void __iter__(void){}
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
