namespace mfem {
%pythonappend ComplexLinearForm::AddDomainIntegrator %{
   self._intg = (lfi_real, lfi_imag)
   if  hasattr(lfi_real, "thisown"): lfi_real.thisown=0
   if  hasattr(lfi_imag, "thisown"): lfi_imag.thisown=0
   lfi_imag.thisown=0
%}  
%pythonappend ComplexLinearForm::AddBoundaryIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[1].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}
%pythonappend ComplexLinearForm::AddBdrFaceIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[1].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend SesquilinearForm::AddDomainIntegrator %{
   self._intg = (bfi_real, bfi_imag)
   if hasattr(bfi_real, "thisown"): bfi_real.thisown=0
   if hasattr(bfi_imag, "thisown"): bfi_imag.thisown=0
%}  
%pythonappend SesquilinearForm::AddBoundaryIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[1].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
%pythonappend SesquilinearForm::AddInteriorFaceIntegrator %{
   self._intg = (bfi_real, bfi_imag)
   if hasattr(bfi_real, "thisown"): bfi_real.thisown=0
   if hasattr(bfi_imag, "thisown"): bfi_imag.thisown=0
%}  
%pythonappend SesquilinearForm::AddBdrFaceIntegrator %{
   self._intg = args
   if  hasattr(args[0], "thisown"): args[1].thisown=0
   if  hasattr(args[1], "thisown"): args[1].thisown=0
%}  
}
