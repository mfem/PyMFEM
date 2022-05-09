namespace mfem {
%pythonappend TransposeIntegrator::TransposeIntegrator %{
    if own_bfi_ == 1:  bfi_.thisown = 0
%}

%pythonappend LumpedIntegrator::LumpedIntegrator %{
    if own_bfi_ == 1:  bfi_.thisown = 0
%}

%pythonappend InverseIntegrator::InverseIntegrator %{
    if own_integ == 1:  integ.thisown = 0
%}

%pythonappend SumIntegrator::SumIntegrator %{
    self.own_integs = own_integs
%}

%pythonappend MixedScalarIntegrator::MixedScalarIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorIntegrator::MixedVectorIntegrator %{
    self._coeff = args
%}

%pythonappend MixedScalarVectorIntegrator::MixedScalarVectorIntegrator %{
    self._coeff = args
%}

%pythonappend MixedScalarMassIntegrator::MixedScalarMassIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorProductIntegrator::MixedVectorProductIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarDerivativeIntegrator::MixedScalarDerivativeIntegrator %{
    self._coeff = args
%}

%pythonappend MixedScalarWeakDerivativeIntegrator::MixedScalarWeakDerivativeIntegrator %{
    self._coeff = args
%}

%pythonappend MixedScalarDivergenceIntegrator::MixedScalarDivergenceIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorDivergenceIntegrator::MixedVectorDivergenceIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarWeakGradientIntegrator::MixedScalarWeakGradientIntegrator %{
    self._coeff = args
%}

%pythonappend MixedScalarCurlIntegrator::MixedScalarCurlIntegrator %{
    self._coeff = args
%}

%pythonappend MixedScalarWeakCurlIntegrator::MixedScalarWeakCurlIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorMassIntegrator::MixedVectorMassIntegrator %{
    self._coeff = args
%}

%pythonappend MixedCrossProductIntegrator::MixedCrossProductIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedDotProductIntegrator::MixedDotProductIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedWeakGradDotIntegrator::MixedWeakGradDotIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedWeakDivCrossIntegrator::MixedWeakDivCrossIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedGradGradIntegrator::MixedGradGradIntegrator %{
    self._coeff = args
%}

%pythonappend MixedCrossGradGradIntegrator::MixedCrossGradGradIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedCurlCurlIntegrator::MixedCurlCurlIntegrator %{
    self._coeff = args
%}

%pythonappend MixedCrossCurlCurlIntegrator::MixedCrossCurlCurlIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedCrossCurlGradIntegrator::MixedCrossCurlGradIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedCrossGradCurlIntegrator::MixedCrossGradCurlIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedWeakCurlCrossIntegrator::MixedWeakCurlCrossIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarWeakCurlCrossIntegrator::MixedScalarWeakCurlCrossIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedCrossGradIntegrator::MixedCrossGradIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedCrossCurlIntegrator::MixedCrossCurlIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarCrossCurlIntegrator::MixedScalarCrossCurlIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarCrossGradIntegrator::MixedScalarCrossGradIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarCrossProductIntegrator::MixedScalarCrossProductIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarWeakCrossProductIntegrator::MixedScalarWeakCrossProductIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedDirectionalDerivativeIntegrator::MixedDirectionalDerivativeIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedGradDivIntegrator::MixedGradDivIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedDivGradIntegrator::MixedDivGradIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedScalarWeakDivergenceIntegrator::MixedScalarWeakDivergenceIntegrator %{
    self._coeff = vq
%}

%pythonappend MixedVectorGradientIntegrator::MixedVectorGradientIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorCurlIntegrator::MixedVectorCurlIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorWeakCurlIntegrator::MixedVectorWeakCurlIntegrator %{
    self._coeff = args
%}

%pythonappend MixedVectorWeakDivergenceIntegrator::MixedVectorWeakDivergenceIntegrator %{
    self._coeff = args
%}

%pythonappend GradientIntegrator::GradientIntegrator %{
    self._coeff = args
%}

%pythonappend DiffusionIntegrator::DiffusionIntegrator %{
    self._coeff = args
%}

%pythonappend MassIntegrator::MassIntegrator %{
    self._coeff = args
%}

%pythonappend BoundaryMassIntegrator::BoundaryMassIntegrator %{
    self._coeff = q
%}

%pythonappend ConvectionIntegrator::ConvectionIntegrator %{
    self._coeff = q
%}

%pythonappend ConservativeConvectionIntegrator::ConservativeConvectionIntegrator %{
    self._coeff = q
%}

%pythonappend GroupConvectionIntegrator::GroupConvectionIntegrator %{
    self._coeff = q
%}

%pythonappend VectorMassIntegrator::VectorMassIntegrator %{
    self._coeff = args
%}

%pythonappend VectorFEDivergenceIntegrator::VectorFEDivergenceIntegrator %{
    self._coeff = args
%}

%pythonappend VectorFEWeakDivergenceIntegrator::VectorFEWeakDivergenceIntegrator %{
    self._coeff = args
%}

%pythonappend VectorFECurlIntegrator::VectorFECurlIntegrator %{
    self._coeff = args
%}

%pythonappend DerivativeIntegrator::DerivativeIntegrator %{
    self._coeff = q
%}

%pythonappend CurlCurlIntegrator::CurlCurlIntegrator %{
    self._coeff = args
%}

%pythonappend VectorCurlCurlIntegrator::VectorCurlCurlIntegrator %{
    self._coeff = args
%}

%pythonappend VectorFEMassIntegrator::VectorFEMassIntegrator %{
    self._coeff = args
%}

%pythonappend VectorDivergenceIntegrator::VectorDivergenceIntegrator %{
    self._coeff = args
%}

%pythonappend DivDivIntegrator::DivDivIntegrator %{
    self._coeff = args
%}

%pythonappend VectorDiffusionIntegrator::VectorDiffusionIntegrator %{
    self._coeff = args
%}

%pythonappend ElasticityIntegrator::ElasticityIntegrator %{
    self._coeff = args
%}

%pythonappend DGTraceIntegrator::DGTraceIntegrator %{
    self._coeff = args
%}

%pythonappend NonconservativeDGTraceIntegrator::NonconservativeDGTraceIntegrator %{
    self._coeff = args
%}

%pythonappend DGDiffusionIntegrator::DGDiffusionIntegrator %{
    self._coeff = args
%}

%pythonappend DGDiffusionBR2Integrator::DGDiffusionBR2Integrator %{
    self._coeff = args
%}

%pythonappend DGElasticityIntegrator::DGElasticityIntegrator %{
    self._coeff = args
%}

%pythonappend ScalarProductInterpolator::ScalarProductInterpolator %{
    self._coeff = sc
%}

%pythonappend ScalarVectorProductInterpolator::ScalarVectorProductInterpolator %{
    self._coeff = sc
%}

%pythonappend VectorScalarProductInterpolator::VectorScalarProductInterpolator %{
    self._coeff = vc
%}

%pythonappend ScalarCrossProductInterpolator::ScalarCrossProductInterpolator %{
    self._coeff = vc
%}

%pythonappend VectorCrossProductInterpolator::VectorCrossProductInterpolator %{
    self._coeff = vc
%}

%pythonappend VectorInnerProductInterpolator::VectorInnerProductInterpolator %{
    self._coeff = vc
%}

%pythonappend SumIntegrator::AddIntegrator %{
   if self.own_integs == 1: integ.thisown = 0
%}
}
