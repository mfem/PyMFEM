// Used in array.i
//  INTPTR_SIZE_IN(int *data_, int asize)
//  DOUBLEPTR_SIZE_IN(double *data_, int asize)
//
//  This typemap accept either list or [pointer *, size]
%define XXXPTR_SIZE_IN(P1, P2, XXX)  
%typemap(in) (P1, P2) (XXX * temp_ptr,
		       //bool ptr_given=false,
			bool is_tuple=false,
			PyObject *s1,
			PyObject *s2){
  if (!PyList_Check($input)) {
      if (!PyTuple_Check($input)) {
         PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
         return NULL;
      } else {
	is_tuple = true;
      }
  }
  s1 = (is_tuple) ? PyTuple_GetItem($input, 0) : PyList_GetItem($input,0);
  s2 = (is_tuple) ? PyTuple_GetItem($input, 1) : PyList_GetItem($input,1);    
  SWIG_ConvertPtr(s1, (void **) &temp_ptr,
		  $descriptor(const XXX *), 0 |0);
  $1 = temp_ptr;
  $2 = PyLong_AsLong(s2);    
}

%typemap(typecheck) (P1, P2) {
  XXX * temp_ptr;
  $1 = 0;
  if (PyList_Check($input)){
     int size = PyList_Size($input);    
     if (size == 2){
          PyObject *s1 = PyList_GetItem($input,0);       
          if (SWIG_ConvertPtr(s1, (void **) &temp_ptr,
   			 $descriptor(const XXX *), 0 |0) == -1) {
 	    //
         } else {
            $1 = 1;
         }
     }
  }
  if (PyTuple_Check($input)){
     int size = PyTuple_Size($input);        
     if (size == 2){
          PyObject *s1 = PyTuple_GetItem($input,0);       
          if (SWIG_ConvertPtr(s1, (void **) &temp_ptr,
   			 $descriptor(const XXX *), 0 |0) == -1) {
 	    //
         } else {
            $1 = 1;
         }
     }
  }
}

%typemap(newfree) (P1, P2) {
   if ($1) free($1);
}
%enddef
   /*
%define DOUBLEPTR_SIZE_IN(P1, P2)
%typemap(in) (P1, P2) (double * temp_ptr,
			bool ptr_given=false,
			bool is_tuple=false,
			PyObject *s1,
			PyObject *s2){
  int i;
  if (!PyList_Check($input)) {
      if (!PyTuple_Check($input)) {
         PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
         return NULL;
      } else {
	is_tuple = true;
      }
  }
  $2 = (is_tuple) ? PyTuple_Size($input) : PyList_Size($input);  
  if ($2 == 2){
    s1 = (is_tuple) ? PyTuple_GetItem($input, 0) : PyList_GetItem($input,0);
    s2 = (is_tuple) ? PyTuple_GetItem($input, 1) : PyList_GetItem($input,1);    
    if (SWIG_ConvertPtr(s1, (void **) &temp_ptr,
			$descriptor(const double *), 0 |0) == -1) {
      ptr_given=false;
    } else {
      ptr_given=true;
    }
  }
  if (! ptr_given){
    $1 = (double *) malloc(($2)*sizeof(double));
    for (i = 0; i < $2; i++) {
      PyObject *s = (is_tuple) ? PyTuple_GetItem($input, i) : PyList_GetItem($input,i);      
      if (PyInt_Check(s)) {
	$1[i] = (double)PyInt_AsLong(s);
      } else if (PyFloat_Check(s)) {	
        $1[i] = (double)PyFloat_AsDouble(s);	
      } else {
	free($1);
	PyErr_SetString(PyExc_ValueError, "List items must be float");
	return NULL;
      }
    }
  } else {
    $1 = temp_ptr;
    $2 = PyLong_AsLong(s2);    
  }
}
%typemap(typecheck) (P1, P2) {
  $1 = 0;
  if (PyList_Check($input)){
    $1 = 1;
  }
  if (PyTuple_Check($input)){
    $1 = 1;
  }
}

%typemap(newfree) (P1, P2) {
   if ($1) free($1);
}
%enddef

   */
