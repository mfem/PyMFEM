namespace mfem {
  // serial
%pythonappend ComplexLinearForm::AddDomainIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend ComplexLinearForm::AddBoundaryIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}
%pythonappend ComplexLinearForm::AddBdrFaceIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend SesquilinearForm::AddDomainIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend SesquilinearForm::AddBoundaryIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend SesquilinearForm::AddInteriorFaceIntegrator %{
   self._intg = (bfi_real, bfi_imag)
   if hasattr(bfi_real, "thisown"): bfi_real.thisown=0
   if hasattr(bfi_imag, "thisown"): bfi_imag.thisown=0
%}  
%pythonappend SesquilinearForm::AddBdrFaceIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}
  // parallel  
%pythonappend ParComplexLinearForm::AddDomainIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend ParComplexLinearForm::AddBoundaryIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}
%pythonappend ParComplexLinearForm::AddBdrFaceIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend ParSesquilinearForm::AddDomainIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend ParSesquilinearForm::AddBoundaryIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend ParSesquilinearForm::AddInteriorFaceIntegrator %{
   self._intg = (bfi_real, bfi_imag)
   if hasattr(bfi_real, "thisown"): bfi_real.thisown=0
   if hasattr(bfi_imag, "thisown"): bfi_imag.thisown=0
%}  
%pythonappend ParSesquilinearForm::AddBdrFaceIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[0].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
  
}
