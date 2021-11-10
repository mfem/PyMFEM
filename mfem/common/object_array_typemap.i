//
//   Array< T > &  ->  Python list of T
// 
//  ObjectArrayInput(Solver*)
%define ObjectArrayInput(T)
%typemap(in) mfem::Array < ## T ## >& ( T temp_ptr,  mfem::Array< ## T ## > tmp_array) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int l = PyList_Size($input);
  for (i = 0; i < l; i++) {
     PyObject *s1 = PyList_GetItem($input,i);
     if (SWIG_ConvertPtr(s1, (void **) &temp_ptr, $descriptor(const T ), 0 |0) == -1) {
          PyErr_SetString(PyExc_ValueError, "Expecting a list of  T ");       
          return NULL;
     }
  }
  tmp_array.SetSize(l);
  for (i = 0; i < l; i++) {
     PyObject *s1 = PyList_GetItem($input,i);
     SWIG_ConvertPtr(s1, (void **) &temp_ptr, $descriptor(const T ), 0 |0);
     tmp_array[i] = temp_ptr;
  }
  $1 = &tmp_array;
}
%typemap(typecheck) mfem::Array < ## T ## >&  {
  $1 = 0;
  if (PyList_Check($input)){
    $1 = 1;
  }
}
%enddef

//
//   Array<bool> &  ->  Python list of bool
// 
%define BoolArrayInput(T)
%typemap(in) mfem::Array <bool>& ( bool temp_ptr,  mfem::Array<bool> tmp_array){
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int l = PyList_Size($input);
  for (i = 0; i < l; i++) {
     PyObject *s1 = PyList_GetItem($input,i);
     if (! PyBool_Check(s1)){
          PyErr_SetString(PyExc_ValueError, "Expecting a list of bool");       
          return NULL;
     }
  }
  tmp_array.SetSize(l);
  for (i = 0; i < l; i++) {
     PyObject *s1 = PyList_GetItem($input,i);
     int isTrue = PyObject_IsTrue(s1);
     if (isTrue){
       tmp_array[i] = true;
     } else {
       tmp_array[i] = false;
     }
  }
  $1 = &tmp_array;
}
%typemap(typecheck) mfem::Array <bool>&  {
  $1 = 0;
  if (PyList_Check($input)){
    $1 = 1;
  }
}
%enddef
