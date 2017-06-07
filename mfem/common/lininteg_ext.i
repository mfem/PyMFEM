namespace mfem {
%pythonappend BoundaryLFIntegrator::BoundaryLFIntegrator %{
    self._coeff = QG
%}
%pythonappend DomainLFIntegrator::DomainLFIntegrator %{
    self._coeff = args[0]
%}
%pythonappend VectorFEBoundaryTangentLFIntegrator::VectorFEBoundaryTangentLFIntegrator %{
    self._coeff = QG
%}
%pythonappend VectorFEDomainLFIntegrator::VectorFEDomainLFIntegrator %{
   self._coeff = F
%}
}

