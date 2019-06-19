%module(package="mfem._par") table
%{
#include "iostream_typemap.hpp"    
#include "general/table.hpp"
#include <iostream>
#include "numpy/arrayobject.h"    
%}

// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%include "exception.i"
%import "array.i"
%import "ostream_typemap.i"
%import "../common/ignore_common_functions.i"
%import "../common/numpy_int_typemap.i"
%import "../common/exception.i"

%import "mem_manager.i"

%include "general/table.hpp"

%extend mfem::Table {
PyObject* GetRowList(int i) const{
     const int *row = self->GetRow(i);
     int L = self->RowSize(i);
     int kk = 0;
     PyObject *o;
     o = PyList_New(L);
     for (kk =0; kk < L; kk ++){
       PyList_SetItem(o, kk, PyLong_FromLong((long) *(row + kk)));
     }
     return o;
  }
};
