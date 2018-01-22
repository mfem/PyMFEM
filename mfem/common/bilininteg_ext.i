namespace mfem {
%pythonprepend TransposeIntegrator::TransposeIntegrator %{
    if _own_bfi == 1:  _bfi.thisown = 0
%}
%pythonprepend InverseIntegrator::InverseIntegrator %{
    if own_integ == 1:  integ.thisown = 0
%}
%pythonprepend SumIntegrator::AddIntegrator %{
    integ.thisown = 0
%}
%pythonappend MassIntegrator::MassIntegrator %{
    if len(args) > 0: self._coeff = args[0]
%}
%pythonappend DiffusionIntegrator::DiffusionIntegrator %{
    if len(args) > 0: self._coeff = args[0]
%}
%pythonappend CurlCurlIntegrator::CurlCurlIntegrator %{
    if len(args) > 0: self._coeff = args[0]
%}
%pythonappend VectorFEMassIntegrator::VectorFEMassIntegrator %{
    if len(args) > 0: self._coeff = args[0]
%}
%pythonappend MixedVectorGradientIntegrator::MixedVectorGradientIntegrator%{
    if len(args) > 0: self._coeff = args[0]
%}
%pythonappend MixedVectorWeakDivergenceIntegrator::MixedVectorWeakDivergenceIntegrator%{
    if len(args) > 0: self._coeff = args[0]
%}
%pythonappend  MixedDotProductIntegrator::MixedDotProductIntegrator%{
    self._coeff = vq
%}
%pythonappend  MixedWeakGradDotIntegrator::MixedWeakGradDotIntegrator%{
    self._coeff = vq
%}
}  
  
