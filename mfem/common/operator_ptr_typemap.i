//
//   allows to pass OperatorPtr to operator & 
//
%typemap(in, noblock=1) mfem::Operator & (void *argp = 0, int res = 0) {
  // default conversion
  res = SWIG_ConvertPtr($input, &argp, $descriptor(mfem::Operator &), $disown |  0 );
  if (!SWIG_IsOK(res)) {
      // try to treat it as OperatorPtr
      res = SWIG_ConvertPtr($input, &argp, $descriptor(mfem::OperatorHandle *), $disown |  0 );
      if (!SWIG_IsOK(res)) {          
         SWIG_exception_fail(SWIG_ArgError(res), "in method '" "$symname" "', argument "
                       "$argnum"" of type '" "$type""'");
      } else {
	$1 = (mfem::Operator *)(((mfem::OperatorHandle const *)argp)->Ptr());
      }
  }else{
    $1 = (mfem::Operator *)(argp);    
  }
}  

