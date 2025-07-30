%import "element.i"
%import "array.i"
%import "vector.i"
%import "blockvector.i"
%import "bilininteg.i"
%import "lininteg.i"
%import "fespace.i"
%import "blockmatrix.i"
%import "operators.i"
%import "blockmatrix.i"
%import "blockoperator.i"
%import "complex_densemat.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)

%ignore mfem::BlockStaticCondensation::ConvertListToReducedTrueDofs;
%ignore mfem::ComplexBlockStaticCondensation::ConvertListToReducedTrueDofs;

namespace mfem {
%pythonprepend DPGWeakForm::AddTrialIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    di = bfi
    self._integrators.append(di)
    di.thisown=0	     ## ownership is C++ side
   %}
%pythonprepend DPGWeakForm::AddTestIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    di = bfi
    self._integrators.append(di)
    di.thisown=0	     ## ownership is C++ side
   %}
%pythonprepend DPGWeakForm::AddDomainLFIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    di = lfi
    self._integrators.append(di)
    di.thisown=0	     ## ownership is C++ side
   %}
%pythonprepend ComplexDPGWeakForm::AddTrialIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    di = bfi_r
    if di is not None:
        self._integrators.append(di)
        di.thisown=0	     ## ownership is C++ side
    di = bfi_i
    if di is not None:
        self._integrators.append(di)
        di.thisown=0	     ## ownership is C++ side
   %}
%pythonprepend ComplexDPGWeakForm::AddTestIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    di = bfi_r
    if di is not None:
        self._integrators.append(di)
        di.thisown=0	     ## ownership is C++ side
    di = bfi_i
    if di is not None:
        self._integrators.append(di)
        di.thisown=0	     ## ownership is C++ side
   %}
%pythonprepend ComplexDPGWeakForm::AddDomainLFIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    di = lfi_r
    if di is not None:
        self._integrators.append(di)
        di.thisown=0	     ## ownership is C++ side
    di = lfi_i
    if di is not None:
        self._integrators.append(di)
        di.thisown=0	     ## ownership is C++ side
   %}
}

//%constant mfem::real_t detJ_r_function(const mfem::Vector &, mfem::CartesianPML *);

%constant mfem::real_t detJ_r_function(const mfem::Vector &, mfem::CartesianPML *);
%constant mfem::real_t detJ_i_function(const mfem::Vector &, mfem::CartesianPML *);
%constant mfem::real_t abs_detJ_2_function(const mfem::Vector &, mfem::CartesianPML *);
%constant void Jt_J_detJinv_r_function(const mfem::Vector & , mfem::CartesianPML *,
					     mfem::DenseMatrix & );
%constant void Jt_J_detJinv_i_function(const mfem::Vector & , mfem::CartesianPML *,
					     mfem::DenseMatrix &);
%constant void abs_Jt_J_detJinv_2_function(const mfem::Vector &, mfem::CartesianPML *,
						 mfem::DenseMatrix &);
%constant void detJ_Jt_J_inv_r_function(const mfem::Vector &, mfem::CartesianPML *,
					      mfem::DenseMatrix &);
%constant void detJ_Jt_J_inv_i_function(const mfem::Vector &, mfem::CartesianPML *,
					      mfem::DenseMatrix &);
%constant void abs_detJ_Jt_J_inv_2_function(const mfem::Vector &, mfem::CartesianPML *,
						  mfem::DenseMatrix &);


%rename(detJ_r_f) mfem::detJ_r_function;
%rename(detJ_i_f) mfem::detJ_i_function;
%rename(abs_detJ_2_f) mfem::abs_detJ_2_function;
%rename(Jt_J_detJinv_r_f) mfem::Jt_J_detJinv_r_function;
%rename(Jt_J_detJinv_i_f) mfem::Jt_J_detJinv_i_function;
%rename(abs_Jt_J_detJinv_2_f) mfem::abs_Jt_J_detJinv_2_function;
%rename(detJ_Jt_J_inv_r_f) mfem::detJ_Jt_J_inv_r_function;
%rename(detJ_Jt_J_inv_i_f) mfem::detJ_Jt_J_inv_i_function;
%rename(abs_detJ_Jt_J_inv_2_f) mfem::abs_detJ_Jt_J_inv_2_function;


// % handle Array input
%include "../common/typemap_macros.i"
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementSpace*>&,
			 mfem::FiniteElementSpace*)
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::ParFiniteElementSpace*>&,
			 mfem::ParFiniteElementSpace*)
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementCollection*>&,
			 mfem::FiniteElementCollection*)


//DPGWeakForm
%pythonprepend mfem::DPGWeakForm::DPGWeakForm %{
  if len(args) > 0:
     fes_, fecol_ = args
     self._fes = fes_
     self._fecol = fecol_
%}

%pythonprepend mfem::DPGWeakForm::SetSpaces %{
  self._fes = fes_
  self._fecol = fecol_
%}

//ComplexDPGWeakForm
%pythonprepend mfem::ComplexDPGWeakForm::ComplexDPGWeakForm %{
  if len(args) > 0:
     fes_, fecol_ = args
     self._fes = fes_
     self._fecol = fecol_
%}

%pythonprepend mfem::ComplexDPGWeakForm::SetSpaces %{
  self._fes = fes_
  self._fecol = fecol_
%}

//ParDPGWeakForm
%pythonprepend mfem::ParDPGWeakForm::ParDPGWeakForm %{
    if len(args) > 0:
     trial_pfes_, fecol_ = args
     self._fes = trial_pfes_
     self._fecol = fecol_

%}

%pythonprepend mfem::ParDPGWeakForm::SetParSpaces %{
  self._pfes = trial_pfes_
  self._fecol = fecol_
%}

//ParComplexDPGWeakForm
%pythonprepend mfem::ParComplexDPGWeakForm::ParComplexDPGWeakForm %{
    if len(args) > 0:
     trial_pfes_, fecol_ = args
     self._fes = trial_pfes_
     self._fecol = fecol_

%}

%pythonprepend mfem::ParComplexDPGWeakForm::SetParSpaces %{
  self._pfes = trial_pfes_
  self._fecol = fecol_
%}
