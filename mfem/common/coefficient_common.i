namespace mfem {
%pythonappend PWCoefficient::PWCoefficient  %{
    if len(args) > 1:
       self._link = args[1]
%}
%pythonappend GridFunctionCoefficient::GridFunctionCoefficient  %{
    if len(args) > 0:
       self._link = args[0]
%}
%pythonprepend GridFunctionCoefficient::SetGridFunction %{
    self._link = gf
%}
%pythonappend TransformedCoefficient::TransformedCoefficient  %{
    self._link = [x for x in args if isinstance(x, Coefficient)]
%}
%pythonprepend DeltaCoefficient::SetWeight %{
    self._linkw = w
%}
%pythonappend RestrictedCoefficient::RestrictedCoefficient %{
    self._ref_to_c = c_
%}
%pythonappend PWVectorCoefficient::PWVectorCoefficient %{
    if len(args) > 2:
       self._link = args[2]      
%}
%pythonprepend VectorArrayCoefficient::Set %{ 
    c.thisown=0 
%}
%pythonprepend VectorGridFunctionCoefficient::VectorGridFunctionCoefficient %{
    if len(args) > 0:
       self._link = args[0]
%}
%pythonprepend VectorGridFunctionCoefficient::SetGridFunction %{
    self._link = gf
%}
%pythonprepend GradientGridFunctionCoefficient::GradientGridFunctionCoefficient %{
    self._link = gf    
%}
%pythonprepend GradientGridFunctionCoefficient::SetGridFunction %{
    self._link = gf
%}
%pythonprepend CurlGridFunctionCoefficient::CurlGridFunctionCoefficient %{
    self._link = gf        
%}
%pythonprepend CurlGridFunctionCoefficient::SetGridFunction %{
    self._link = gf
%}
%pythonprepend DivergenceGridFunctionCoefficient::DivergenceGridFunctionCoefficient %{
    self._link = gf    
%}
%pythonprepend DivergenceGridFunctionCoefficient::SetGridFunction %{
    self._link = gf
%}
%pythonappend VectorRestrictedCoefficient::VectorRestrictedCoefficient %{
    self._ref_to_vc = vc
%}

%pythonprepend MatrixArrayCoefficient::Set %{ 
    c.thisown=0 
%}
  
%pythonappend MatrixRestrictedCoefficient::MatrixRestrictedCoefficient %{
    self._ref_to_mc = mc
%}

%pythonappend SumCoefficient::SumCoefficient %{
    if len(args) > 0 and isinstance(agrs[0], Coefficient):
          self.linkA = args[0]
    if len(args) > 1 and isinstance(agrs[1], Coefficient):
          self.linkB = args[0]

%}
%pythonappend SumCoefficient::SetACoef %{
    self.linkA = A    
%}    
%pythonappend SumCoefficient::SetBCoef %{
    self.linkB = B    
%}    
  
%pythonappend ProductCoefficient::ProductCoefficient %{
    if len(args) > 0 and isinstance(args[0], Coefficient):
       self._linkA = args[0]
    if len(args) > 1 and isinstance(args[1], Coefficient):
       self._linkB = args[1]
%}  
%pythonappend ProductCoefficient::SetACoef %{
       self._linkA = A    
%}
%pythonappend ProductCoefficient::SetBCoef %{
       self._linkB = B    
%}

  
%pythonappend RatioCoefficient::RatioCoefficient %{
    if len(args) > 0 and isinstance(args[0], Coefficient):
       self._linkA = args[0]
    if len(args) > 1 and isinstance(args[1], Coefficient):
       self._linkB = args[1]
%}  
%pythonappend RatioCoefficient::SetACoef %{
       self._linkA = A    
%}
%pythonappend RatioCoefficient::SetBCoef %{
       self._linkB = B    
%}

%pythonappend PowerCoefficient::PowerCoefficient %{
    self._linkA = A
%}  
%pythonappend PowerCoefficient::SetACoef %{
    self._linkA = A
%}    
  
%pythonappend InnerProductCoefficient::InnerProductCoefficient %{
    self._linkA = A
    self._linkB = B
%}  
%pythonappend InnerProductCoefficient::SetACoef %{
    self._linkA = A
%}    
%pythonappend InnerProductCoefficient::SetBCoef %{
    self._linkB = B
%}    
  

%pythonappend VectorRotProductCoefficient::VectorRotProductCoefficient %{
    self._linkA = A
    self._linkB = B
%}  
%pythonappend VectorRotProductCoefficient::SetACoef %{
    self._linkA = A
%}    
%pythonappend VectorRotProductCoefficient::SetBCoef %{
    self._linkB = B
%}    
  
  
%pythonappend DeterminantCoefficient::DeterminantCoefficient %{
    self.linkA = A
%}    
%pythonappend DeterminantCoefficient::SetACoef%{
    self.linkA = A
%}    

%pythonappend VectorSumCoefficient::VectorSumCoefficient %{
    if len(args) > 1:
        self.linkA = args[0]
        self.linkB = args[1]
    if len(args) > 2:
       if isinstance(agrs[2], Coefficient):
          self.linkAlphaCoef = args[2]
    if len(args) > 3:	    
       if isinstance(agrs[3], Coefficient):
          self.linkBetaCoef = args[3]	
%}
%pythonappend VectorSumCoefficient::SetACoef %{
    self.linkA = A_    
%}    
%pythonappend VectorSumCoefficient::SetBCoef %{
    self.linkB = B_    
%}    
%pythonappend VectorSumCoefficient::SetAlphaCoef %{
    self.linkAlpha = A_    
%}    
%pythonappend VectorSumCoefficient::SetBetaCoef %{
    self.linkBeta = B_    
%}
  
%pythonappend ScalarVectorProductCoefficient::ScalarVectorProductCoefficient %{
    if isinstance(args[0], Coefficient):
        self._linkA = args[0]
    if isinstance(args[1], VectorCoefficient):    
        self._linkB = args[1]
%}    
%pythonappend ScalarVectorProductCoefficient::SetACoef %{
    self._linkA = A
%}    
%pythonappend ScalarVectorProductCoefficient::SetBCoef %{
    self._linkB = B
%}
  
%pythonappend NormalizedVectorCoefficient::NormalizedVectorCoefficient%{
    self._linkA = A
%}    
%pythonappend NormalizedVectorCoefficient::SetACoef %{
    self._linkA = A
%}
  
%pythonappend VectorCrossProductCoefficient::VectorCrossProductCoefficient %{
    self._linkA = A
    self._linkB = B
%}  
%pythonappend VectorCrossProductCoefficient::SetACoef %{
    self._linkA = A
%}    
%pythonappend VectorCrossProductCoefficient::SetBCoef %{
    self._linkB = B
%}    

%pythonappend MatrixVectorProductCoefficient::MatrixVectorProductCoefficient %{
    self._linkA = A
    self._linkB = B
%}    
%pythonappend MatrixVectorProductCoefficient::SetACoef %{
    self._linkA = A
%}    
%pythonappend MatrixVectorProductCoefficient::SetBCoef %{
    self._linkB = B
%}    

%pythonappend MatrixSumCoefficient::MatrixSumCoefficient %{
    self._linkA = A
    self._linkB = B
%}
%pythonappend MatrixSumCoefficient::SetACoef %{
    self._linkA = A
%}
%pythonappend MatrixSumCoefficient::SetBCoef %{
    self._linkB = B
%}
%pythonappend MatrixProductCoefficient::MatrixProductCoefficient %{
    self._linkA = A
    self._linkB = B
%}    
%pythonappend MatrixProductCoefficient::SetACoef %{
    self._linkA = A
%}    
%pythonappend MatrixProductCoefficient::SetBCoef %{
    self._linkB = B
%}    
%pythonappend ScalarMatrixProductCoefficient::ScalarMatrixProductCoefficient %{
    if isinstance(args[0], Coefficient):
        self._linkA = args[0]
    if isinstance(args[1], MatrixCoefficient):    
        self._linkB = args[1]
%}
%pythonappend ScalarMatrixProductCoefficient::SetACoef %{
    self._linkA = A
%}
%pythonappend ScalarMatrixProductCoefficient::SetBCoef %{
    self._linkB = B
%}
									      
%pythonappend TransposeMatrixCoefficient::TransposeMatrixCoefficient%{
   self._link = A
%}
%pythonappend TransposeMatrixCoefficient::SetACoef %{
   self._link = A  
%}
%pythonappend InverseMatrixCoefficient::InverseMatrixCoefficient%{
   self._link = A
%}
%pythonappend InverseMatrixCoefficient::SetACoef%{
   self._link = A
%}
%pythonappend OuterProductCoefficient::OuterProductCoefficient%{
   self._linkA = A
   self._linkB = B
%}
%pythonappend OuterProductCoefficient::SetACoef%{
   self._linkA = A
%}
%pythonappend OuterProductCoefficient::SetBCoef%{
   self._linkB = B
%}
  
%pythonappend CrossCrossCoefficient::CrossCrossCoefficient%{
    if isinstance(args[0], Coefficient):
        self._linkA = args[0]
    if isinstance(args[1], VectorCoefficient):    
        self._linkK = args[1]
%}
%pythonappend CrossCrossCoefficient::SetACoef%{
    self._linkA = A
%}
%pythonappend CrossCrossCoefficient::SetKCoef%{
    self._linkK = K
%}

%pythonappend QuadratureFunctionCoefficient::QuadratureFunctionCoefficient %{
    self._linkQF = qf
%}
  
%pythonappend VectorQuadratureFunctionCoefficient::VectorQuadratureFunctionCoefficient %{
    self._linkQF = qf
%}
  

}
