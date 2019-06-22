//  conversion of Int (can handle numpy int)
%typemap(in) int {
  if ((PyArray_PyIntAsInt($input) == -1) && PyErr_Occurred()) {
    SWIG_exception_fail(SWIG_TypeError, "Input must be integer");
  };  
  $1 = PyArray_PyIntAsInt($input);
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) int {
  if ((PyArray_PyIntAsInt($input) == -1) && PyErr_Occurred()) {
    PyErr_Clear();
    $1 = 0;
  } else {
    $1 = 1;    
  }
}

