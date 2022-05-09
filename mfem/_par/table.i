%module(package="mfem._par") table

%feature("autodoc", "1");

%{
#include <fstream>  
#include <iostream>
#include "mfem.hpp"  
#include "numpy/arrayobject.h"      
#include "../common/io_stream.hpp"
#include "general/table.hpp"
%}

// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%include "exception.i"
%import "array.i"
%import "../common/ignore_common_functions.i"
%import "../common/numpy_int_typemap.i"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

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
/*
  void Print(std::ostream & out = mfem::out, int width = 4) const;
  void PrintMatlab(std::ostream & out) const;
  void Save(std::ostream &out) const;
*/
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(Table, Print)
OSTREAM_ADD_DEFAULT_FILE(Table, PrintMatlab)
OSTREAM_ADD_DEFAULT_STDOUT_FILE(Table, Save)
#endif
