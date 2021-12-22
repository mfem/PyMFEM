/*

   complex_operator.i

*/
%module(package="mfem._par") complex_operator
%feature("autodoc", "1");
%{
#include "linalg/complex_operator.hpp"
#include "numpy/arrayobject.h"
#include "linalg/hypre.hpp"  
#include "pyoperator.hpp"     
  %}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
%import "hypre.i"
#endif

%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "operators.i"
%import "sparsemat.i"

namespace mfem { 
%pythonprepend ComplexOperator::ComplexOperator %{
    self._parts = [Op_Real, Op_Imag]
    if ownReal:
       assert Op_Real.thisown != 0, "Real Op is not owned by passed object"
       Op_Real.thisown = 0
    if ownImag:
       assert Op_Real.thisown != 0, "Imag Op is not owned by passed object"     
       Op_Imag.thisown = 0
    %}
%ignore ComplexOperator::ComplexOperator(Operator * Op_Real, Operator * Op_Imag,
                        bool ownReal, bool ownImag,
                        Convention convention = HERMITIAN);
} //end of namespace

%include "linalg/complex_operator.hpp"

%extend mfem::ComplexOperator{
   ComplexOperator(mfem::Operator *Op_Real,
		   mfem::Operator *Op_Imag,
		   bool ownReal = false,
		   bool ownImag = false,
                   bool hermitan  = true){
   if (hermitan){
     return  new mfem::ComplexOperator(Op_Real, Op_Imag, ownReal, ownImag,
				       mfem::ComplexOperator::HERMITIAN);
   } else {
     return  new mfem::ComplexOperator(Op_Real, Op_Imag, ownReal, ownImag,
				       mfem::ComplexOperator::BLOCK_SYMMETRIC);
   }
}
};
