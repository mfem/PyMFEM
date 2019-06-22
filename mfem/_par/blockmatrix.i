%module(package="mfem._par") blockmatrix

%{
#include <fstream>  
#include <iostream>
#include "linalg/blockmatrix.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "iostream_typemap.hpp"    
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "matrix.i"
%import "sparsemat.i"
%import "ostream_typemap.i"
%import "../common/ignore_common_functions.i"

%pythonappend mfem::BlockMatrix::BlockMatrix %{
from mfem.par import intArray  
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

%include "linalg/blockmatrix.hpp"
