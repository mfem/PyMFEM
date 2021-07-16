// integer array output with known length (tuple)
%define INTARRAY_OUT_TO_TUPLE(type_name, l)
%typemap(out) type_name{
  $result = PyTuple_New(l);
  for(int i = 0; i < l; i++) {
    PyTuple_SetItem($result, i, PyInt_FromLong($1[i]));
  }
}
%enddef

// integer array output with known length (list)
%define INTARRAY_OUT_TO_LIST(type_name, l)
%typemap(out) type_name{
  $result = PyList_New(l);
  for(int i = 0; i < l; i++) {
    PyList_SetItem($result, i, PyInt_FromLong($1[i]));
  }
}
%enddef

// integer output as int point
%define INTARRAY_OUT_TO_INT(type_name)
%typemap(out) type_name{
  $result = PyInt_FromLong($1[0]);
}
%enddef

// wrap integer with  default -1
%define INT_DEFAULT_NEGATIVE_ONE(type_name)
%typemap(in) (type_name) {
  if (PyInt_Check($input)) {
     $1 = PyInt_AsLong($input);
  } else if ((PyArray_PyIntAsInt($input) != -1) || !PyErr_Occurred()) {
     $1 = PyArray_PyIntAsInt($input);
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a integer");
    return NULL;
  }
}
%typemap(typecheck) (type_name) {
  if (PyInt_Check($input)) {
    $1 = 1;
  } else if ((PyArray_PyIntAsInt($input) != -1) || !PyErr_Occurred()) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}
%enddef

// known length of list/tuple to int pointer
%define LIST_TO_INTARRAY_IN(type_name, l)
%typemap(in) type_name (int temp[l]){
  if (PyList_Check($input)) {
     int ll = PyList_Size($input);
     for (int i = 0; i < ll; i++) {
        PyObject *s = PyList_GetItem($input,i);
        temp[i] = (int)PyInt_AsLong(s);
     }
  } else if (PyTuple_Check($input)) {
     int ll = PyTuple_Size($input);
     for (int i = 0; i < ll; i++) {
        PyObject *s = PyTuple_GetItem($input,i);
        temp[i] = (int)PyInt_AsLong(s);
     }
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
    return NULL;
  }

  $1 = temp;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) (type_name) {
  $1 = 0;
  if (PyList_Check($input)){
    if  (PyList_Size($input) == 2){
        $1 = 1;
    }
  }
  if (PyTuple_Check($input)){
    if  (PyTuple_Size($input) == 2){
        $1 = 1;
    }
  }
  
}
%enddef

// int pointer input for single int
%define INT_TO_INTARRAY_IN(type_name)
%typemap(in) type_name (int temp){
  if (PyInt_Check($input)) {
     temp = PyInt_AsLong($input);
  } else if ((PyArray_PyIntAsInt($input) != -1) || !PyErr_Occurred()) {
     temp = PyArray_PyIntAsInt($input);
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a integer");
    return NULL;
  }
  $1 = &temp;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) (type_name) {
  if (PyInt_Check($input)) {
    $1 = 1;
  } else if ((PyArray_PyIntAsInt($input) != -1) || !PyErr_Occurred()) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}
%enddef


// double pointer with undefined length
//  this macro accespt either
//     numpy array: data is managed by numpy
//     or
//     list: list is copyed value will be managed by MFEM
%define ARRAY_TO_DOUBLEARRAY_IN(type_name)
%typemap(in) (type_name){
  int i, si;
  if (SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor, $disown|0) != -1){
	      
  }
  else if (PyArray_Check($input)){
    $1 = (double *) PyArray_DATA((PyArrayObject *)$input);
	 //     $1 = (double *) PyArray_DATA($input);
  }
  else {
     if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_ValueError, "Expecting a list");
        return NULL;
     }
     si = PyList_Size($input);
     $1 = (double *) malloc((si)*sizeof(double));
     for (i = 0; i < si; i++) {
        PyObject *s = PyList_GetItem($input,i);
        if (PyInt_Check(s)) {
            $1[i] = (double)PyFloat_AsDouble(s);
        } else if (PyFloat_Check(s)) {
            $1[i] = (double)PyFloat_AsDouble(s);
        } else {
            PyErr_SetString(PyExc_ValueError, "List items must be integer/float");
            return NULL;
        }
     }
  }

}

%typemap(typecheck ) (type_name){
   if (SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor, 1) != -1){
      $1 = 1;
   }
   else if (PyList_Check($input)){
      $1 = 1;
   }
   else if (PyArray_Check($input)){
      $1 = 1;
   }
   else {
      $1 = 0;
   }
}
%enddef

%define LIST_TO_OBJARRAY_IN(type_name, OBJTYPE)
%typemap(in) type_name (){
  //  List/Tuple -> OBJTYPE
  int res = 0;
  if (PyList_Check($input)) {
     int ll = PyList_Size($input);
     $1 = new mfem::Array<OBJTYPE>(ll);
     for (int i = 0; i < ll; i++) {
       OBJTYPE ttt;
       PyObject *s = PyList_GetItem($input,i);
       res = SWIG_ConvertPtr(s, (void **) &ttt,
			     $descriptor(OBJTYPE),
			     0);
       if (!SWIG_IsOK(res)) { 	 
             return NULL;
       }	
       $1[i] = ttt;
     }
  } else if (PyTuple_Check($input)) {
     int ll = PyTuple_Size($input);
     $1 = new mfem::Array<OBJTYPE>(ll);     
     for (int i = 0; i < ll; i++) {
       OBJTYPE ttt;
       PyObject *s = PyTuple_GetItem($input,i);
       res = SWIG_ConvertPtr(s, (void **) &ttt,
			     $descriptor(OBJTYPE),
			     0);
       if (!SWIG_IsOK(res)) {
	 return NULL;
       }	
       $1[i] = ttt;
     }
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
    return NULL;
  }
     //$1 = temp;
}
%typemap(freearg) type_name{
  if ($1 != 0){
    //delete $1;
  }
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) (type_name) {
  $1 = 0;
  if (PyList_Check($input)){
      $1 = 1;
  }
  if (PyTuple_Check($input)){
     $1 = 1;
  }
}
%enddef

%define LIST_TO_MFEMOBJ_ARRAY_IN(type_name, OBJTYPE)
%typemap(in) type_name (){
  //  List/Tuple -> OBJTYPE
  int res = 0;
  if (PyList_Check($input)) {
     int ll = PyList_Size($input);
     $1 = new mfem::Array<OBJTYPE>(ll);
     for (int i = 0; i < ll; i++) {
       OBJTYPE ttt;
       PyObject *s = PyList_GetItem($input,i);
       if (s == Py_None){
	 ttt = NULL;
       } else {
         res = SWIG_ConvertPtr(s, (void **) &ttt,
			     $descriptor(OBJTYPE),
			     0);
       }
       if (!SWIG_IsOK(res)) { 	 
             return NULL;
       }	
       $1[0][i] = ttt;
     }
  } else if (PyTuple_Check($input)) {
     int ll = PyTuple_Size($input);
     $1 = new mfem::Array<OBJTYPE>(ll);     
     for (int i = 0; i < ll; i++) {
       OBJTYPE ttt;
       PyObject *s = PyTuple_GetItem($input,i);
       if (s == Py_None){
	 ttt = NULL;
       } else {
         res = SWIG_ConvertPtr(s, (void **) &ttt,
			     $descriptor(OBJTYPE),
			     0);
       }
       if (!SWIG_IsOK(res)) {
	 return NULL;
       }	
       $1[0][i] = ttt;
     }
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
    return NULL;
  }
     //$1 = temp;
}
%typemap(freearg) type_name{
  if ($1 != 0){
    delete $1;
  }
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER) (type_name) {
  $1 = 0;
  if (PyList_Check($input)){
     $1 = 1;
  }
  if (PyTuple_Check($input)){
     $1 = 1;
  }
}
%enddef

