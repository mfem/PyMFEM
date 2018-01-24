//  conversion of Int (can handle numpy int)
%typemap(in) int {
  if (PyArray_PyIntAsInt($input) == -1){
    SWIG_exception_fail(SWIG_TypeError, "Input must be integer");
  };  
  $1 = PyInt_AsLong($input);
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) int {
  if (PyArray_PyIntAsInt($input)   != -1){
    $1 = 1;
  } else {
    $1 = 0;
  }
}
