%module(package="mfem._ser") sparsesmoothers
%{
#include "linalg/sparsesmoothers.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "sparsemat.i"
%import "matrix.i"
%import "../common/exception.i"

%include "linalg/sparsesmoothers.hpp"

 // this should not be necessary.
 // dynamic_cast at DSO boundary issue...
/*
%extend mfem::SparseSmoother{
  //  void mfem::SparseSmoother::SetOperator(const SparseMatrix &a) {
  void SetOperator(mfem::SparseMatrix &a) {  
    std::cout << "this is called\n";
    std::cout << std::to_string(a.Height()) << "\n";
    mfem::SparseMatrix *b = new mfem::SparseMatrix(a);
    $self->SetOperator(*b);
  }
};
*/




  
