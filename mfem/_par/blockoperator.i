%module(package="mfem._par") blockoperator

%{
#include "linalg/blockoperator.hpp"
#include "numpy/arrayobject.h"
#include "linalg/operator.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/hypre.hpp"    
#include "pyoperator.hpp"       
%}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "operators.i"

%pythonappend mfem::BlockOperator::BlockOperator %{
from mfem.par import intArray  
if len(args) == 1:
   if isinstance(args[0], intArray):
       self._offsets = args[0]
if len(args) == 2:
   if (isinstance(args[0], intArray) and
       isinstance(args[1], intArray)):
       self._offsets = (args[0], args[1])
%}
%pythonappend mfem::BlockOperator::SetDiagonalBlock %{
  if not hasattr(self, '_linked_op'):
     self._linked_op = {}
  self._linked_op[iblock, iblock] = op
%}  
%pythonappend mfem::BlockOperator::SetBlock %{
  if not hasattr(self, '_linked_op'):
     self._linked_op = {}
  self._linked_op[iRow, iCol] = op
%}

%pythonappend mfem::BlockDiagonalPreconditioner::BlockDiagonalPreconditioner %{
   self._offsets = offsets
%}

%pythonappend mfem::BlockDiagonalPreconditioner::SetDiagonalBlock %{
  if not hasattr(self, '_linked_op'):
     self._linked_op = {}
  self._linked_op[iblock, iblock] = op
%}

%pythonappend mfem::BlockLowerTriangularPreconditioner::BlockLowerTriangularPreconditioner %{
   self._offsets = offsets
%}
%pythonappend mfem::BlockLowerTriangularPreconditioner::SetDiagonalBlock %{
  if not hasattr(self, '_linked_op'):
     self._linked_op = {}
  self._linked_op[iblock, iblock] = op
%}
%pythonprepend mfem::BlockLowerTriangularPreconditioner::SetBlock %{
  if not (iRow > iCol):
      raise ValueError("can not set upper triangle")
%}
%pythonappend mfem::BlockLowerTriangularPreconditioner::SetBlock %{
   if not hasattr(self, '_linked_op'):
      self._linked_op = {}

   self._linked_op[iRow, iCol] = op
%}  

%inline %{
  mfem::BlockOperator *Opr2BlockOpr(mfem::Operator *op) {
    return dynamic_cast<mfem::BlockOperator*>(op);
  }
  mfem::SparseMatrix *Opr2SparseMat(mfem::Operator *op) {
    return dynamic_cast<mfem::SparseMatrix*>(op);
  }
  mfem::HypreParMatrix *Opr2HypreParMat(mfem::Operator *op) {
    return dynamic_cast<mfem::HypreParMatrix*>(op);
  }
  
%}
%include "linalg/blockoperator.hpp"
