%module(package="mfem._ser", directors="1") solvers

%{
#include "linalg/handle.hpp"  
#include "linalg/matrix.hpp"
#include "linalg/sparsemat.hpp"
#include "linalg/solvers.hpp"
#include "pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "numpy/arrayobject.h"
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"

%include "exception.i"
%import "globals.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "sparsemat.i"
%import "../common/exception_director.i"
%import "../common/operator_ptr_typemap.i"
%import "../common/exception_director.i"

%ignore mfem::IterativeSolverMonitor::SetIterativeSolver;
%feature("director") mfem::IterativeSolverMonitor;
%feature("director") mfem::PyIterativeSolver;

%include "linalg/solvers.hpp"
%include "../common/pysolvers.hpp"

#ifdef MFEM_USE_SUITESPARSE
%rename($ignore, %$isfunction) "";
%rename($ignore, %$isclass) "";
%include <umfpack.h>
%include <klu.h>

// access to UMFPackSolver.Control and Info
%pythoncode %{
UMFPackSolver._Control = UMFPackSolver.Control
def getControl(self):
    import ctypes
    mem = (ctypes.c_double * UMFPACK_CONTROL).from_address(int(self._Control))
    pointer = ctypes.pointer(mem)
    return pointer.contents
UMFPackSolver.Control = property(getControl)

UMFPackSolver._Info = UMFPackSolver.Info
def getInfo(self):
    import ctypes
    mem = (ctypes.c_double * UMFPACK_INFO).from_address(int(self._Info))
    pointer = ctypes.pointer(mem)
    return pointer.contents
UMFPackSolver.Info = property(getInfo)
%}
#endif

%inline %{
namespace mfem{
  void PyIterativeSolver::Mult(const Vector &b, Vector &x) const{
    mfem_error("Mult is not implemented");
  }
  void PyIterativeSolver::MultTranspose(const Vector &b, Vector &x) const{
    mfem_error("MultTranspose is not implemented");    
  }
  void PyIterativeSolver::SetPreconditioner(Solver &pr){
    mfem_error("SetPreconditioner is not implemented");      
  }
  /// Also calls SetOperator for the preconditioner if there is one
  void PyIterativeSolver::SetOperator(const Operator &op){
    mfem_error("SetOperator is not implemented");        
  }
} /* end of namespace */
%}
