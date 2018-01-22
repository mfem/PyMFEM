%module table
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

%import array.i
%import ostream_typemap.i
%import "../common/ignore_common_functions.i"

//  conversion of Int (can handle numpy int)
%typemap(in) int {
  PyArray_PyIntAsInt($input);  
  $1 = PyInt_AsLong($input);
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) int {
  if (PyArray_PyIntAsInt($input)   != -1){
    $1 = 1;
  } else {
    $1 = 0;
  }
}

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
