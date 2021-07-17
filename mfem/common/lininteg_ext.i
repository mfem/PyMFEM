namespace mfem {
%pythonappend LinearFormIntegrator::LinearFormIntegrator %{
    self._coeff = args
%}
%pythonappend DeltaLFIntegrator::DeltaLFIntegrator %{
    self._coeff = args
%}
%pythonappend DomainLFIntegrator::DomainLFIntegrator %{
    self._coeff = args
%}
%pythonappend DomainLFGradIntegrator::DomainLFGradIntegrator %{
    self._coeff = QF
%}
%pythonappend BoundaryLFIntegrator::BoundaryLFIntegrator %{
    self._coeff = QG
%}
%pythonappend BoundaryNormalLFIntegrator::BoundaryNormalLFIntegrator %{
    self._coeff = QG
%}
%pythonappend BoundaryTangentialLFIntegrator::BoundaryTangentialLFIntegrator %{
    self._coeff = QG
%}
%pythonappend VectorDomainLFIntegrator::VectorDomainLFIntegrator %{
    self._coeff = QF
%}
%pythonappend VectorBoundaryLFIntegrator::VectorBoundaryLFIntegrator %{
    self._coeff = QG
%}
%pythonappend VectorFEDomainLFIntegrator::VectorFEDomainLFIntegrator %{
    self._coeff = F
%}
%pythonappend VectorFEDomainLFCurlIntegrator::VectorFEDomainLFCurlIntegrator %{
    self._coeff = F
%}
%pythonappend VectorFEDomainLFDivIntegrator::VectorFEDomainLFDivIntegrator %{
    self._coeff = QF
%}
%pythonappend VectorBoundaryFluxLFIntegrator::VectorBoundaryFluxLFIntegrator %{
    self._coeff = (f, ir)
%}
%pythonappend VectorFEBoundaryFluxLFIntegrator::VectorFEBoundaryFluxLFIntegrator %{
    self._coeff = args
%}
%pythonappend VectorFEBoundaryTangentLFIntegrator::VectorFEBoundaryTangentLFIntegrator %{
    self._coeff = QG
%}
%pythonappend BoundaryFlowIntegrator::BoundaryFlowIntegrator %{
    self._coeff = args
%}
%pythonappend DGDirichletLFIntegrator::DGDirichletLFIntegrator %{
    self._coeff = args
%}
%pythonappend DGElasticityDirichletLFIntegrator::DGElasticityDirichletLFIntegrator %{
    self._coeff = uD_
%}
%pythonappend VectorQuadratureLFIntegrator::VectorQuadratureLFIntegrator %{
    self._coeff = args
%}
%pythonappend QuadratureLFIntegrator::QuadratureLFIntegrator %{
    self._coeff = args
%}
}