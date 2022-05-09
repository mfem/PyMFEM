%module(package="mfem._ser", directors="1")  bilinearform
%{
#include "fem/bilinearform.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
using namespace mfem;  
%}

%init %{
import_array();
%}
%include "exception.i"

%import "globals.i"
%import "mem_manager.i"
%import "array.i"
%import "fespace.i"
%import "fe_coll.i"
%import "intrules.i"
%import "matrix.i"
%import "vector.i"
%import "densemat.i"
%import "sparsemat.i"
%import "lininteg.i"
%import "eltrans.i"
%import "bilininteg.i"
%import "linearform.i"
%import "gridfunc.i"
%include "../common/exception_director.i"

%feature("director") mfem::BilinearForm;

 //%include "fem/coefficient.hpp"
namespace mfem { 
%pythonprepend BilinearForm::AddDomainIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    bfi = args[0]
    self._integrators.append(bfi)
    self.UseExternalIntegrators()
    #bfi.thisown=0
    %}
%pythonprepend BilinearForm::AddBoundaryIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    bfi = args[0]
    self._integrators.append(bfi)
    self.UseExternalIntegrators()
    #bfi.thisown=0 
   %} 
%pythonprepend BilinearForm::AddBdrFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    bfi = args[0]
    self._integrators.append(bfi)
    self.UseExternalIntegrators()
    #bfi.thisown=0
   %} 
%pythonprepend BilinearForm::AddInteriorFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    self.UseExternalIntegrators()
    #bfi.thisown=0
   %}
%pythonappend BilinearForm::SpMat %{
    if not hasattr(self, "_spmat"): self._spmat = []
    self._spmat.append(val)
    val.thisown=0
   %}
%pythonappend BilinearForm::EnableHybridization %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(constr_integ)
    # constr_integ is deleted by Hybridization destructor
    constr_integ.thisown = 0
   %} 
%pythonprepend MixedBilinearForm::AddDomainIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0
   %}
%pythonprepend MixedBilinearForm::AddBoundaryIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    bfi = args[0]
    self._integrators.append(bfi)
    bfi.thisown=0
   %} 
%pythonprepend MixedBilinearForm::AddTraceFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0
   %}
%pythonprepend MixedBilinearForm::AddBdrTraceFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    bfi = args[0]
    self._integrators.append(bfi)
    self.UseExternalIntegrators()
    #bfi.thisown=0 
   %} 
%pythonappend MixedBilinearForm::SpMat %{
    if not hasattr(self, "_spmat"): self._spmat = []
    self._spmat.append(val)
    val.thisown=0
   %}
%pythonprepend DiscreteLinearOperator::AddDomainInterpolator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(di)
    di.thisown=0
   %} 
%pythonprepend DiscreteLinearOperator::AddTraceFaceInterpolator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(di)
    di.thisown=0
    %}
} //end of namespace

/* define BilinearFormPtrArray */
%include "../common/deprecation.i"
DEPRECATED_METHOD(mfem::BilinearForm::GetFES())

/* define IntegrationPointArray and IntegrationRulePtrArray */
%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::BilinearForm *, 1)

%import "../common/data_size_typemap.i"
XXXPTR_SIZE_IN(mfem::BilinearForm **data_, int asize, mfem::BilinearForm *)

%import "../common/array_instantiation_macro.i"
IGNORE_ARRAY_METHODS(mfem::BilinearForm *)
INSTANTIATE_ARRAY0(BilinearForm *, BilinearForm, 1)

%include "fem/bilinearform.hpp"

// instatitate template methods 
%define FORM_SYSTEM_MATRIX_WRAP(OsType)
%template(FormLinearSystem) mfem::BilinearForm::FormLinearSystem<mfem:: ## OsType>;
%template(FormSystemMatrix) mfem::BilinearForm::FormSystemMatrix<mfem:: ## OsType>;
%enddef

FORM_SYSTEM_MATRIX_WRAP(SparseMatrix)

 /*
#ifdef MFEM_USE_MPI
  FORM_SYSTEM_MATRIX_WRAP(mfem::HypreParMatrix)
#endif
 */
  
#ifdef MFEM_USE_PETSC
  FORM_SYSTEM_MATRIX_WRAP(mfem::PetscParMatrix)
#endif  
