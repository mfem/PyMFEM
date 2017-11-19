%module nonlinearform
%{
#include "fem/linearform.hpp"    
#include "fem/nonlininteg.hpp"
#include "fem/nonlinearform.hpp"
#include "pyoperator.hpp"     
%}
/*
%init %{
import_array();
%}
*/
%import operators.i
%import fespace.i
%import nonlininteg.i

 //%include "fem/coefficient.hpp"
namespace mfem { 
%pythonprepend NonlinearForm::AddDomainIntegrator %{
#    if not hasattr(self, "_integrators"): self._integrators = []
#    self._integrators.append(nlfi)
    nlfi.thisown=0 
%}
}

%include "fem/nonlinearform.hpp"

%extend mfem::NonlinearForm{
  mfem::SparseMatrix *GetGradientMatrix(const mfem::Vector &x) const
  {
      return  dynamic_cast<mfem::SparseMatrix *>(&$self->GetGradient(x));
  }
};
