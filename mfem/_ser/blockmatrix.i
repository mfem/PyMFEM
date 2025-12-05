//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") blockmatrix

%{
#include <fstream>
#include <iostream>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/io_stream.hpp"
%}
// initialization required to return numpy array from SWIG
%init %{
import_array1(-1);
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "matrix.i"
%warn "Coming here"
%import "sparsemat.i"
%import "../common/exception.i"
%import "../common/ignore_common_functions.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%pythonappend mfem::BlockMatrix::BlockMatrix %{
from mfem.ser import intArray
if len(args) == 1:
   if isinstance(args[0], intArray):
       self._offsets = args[0]
if len(args) == 2:
   if (isinstance(args[0], intArray) and
       isinstance(args[1], intArray)):
       self._offsets = (args[0], args[1])
%}

%pythonappend mfem::BlockMatrix::SetBlock %{
  if not hasattr(self, '_linked_mat'):
     self._linked_mat = {}
  self._linked_mat[i, j] = mat
%}

%inline %{
  mfem::BlockMatrix *Opr2BlockMatrix(mfem::Operator *op) {
    return dynamic_cast<mfem::BlockMatrix*>(op);
  }
%}
%include "linalg/blockmatrix.hpp"

 /*
void PrintMatlab(std::ostream & os = mfem::out) const;
 */
#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(BlockMatrix, PrintMatlab)
#endif

