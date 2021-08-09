// this typemap accept list/const */vector/numpy float array for double *  (used in mesh.i)
// 
//CONST_DOUBLEPTR_IN(const double *)
%define CONST_DOUBLEPTR_IN(TYPENAME)
%typemap(in) (TYPENAME) (mfem::Vector *temp_vec, double * temp_ptr){
  int i;
  if (!PyList_Check($input)) {
    if (SWIG_ConvertPtr($input, (void **) &temp_ptr, $descriptor(const double *), 0 |0) == -1) {
       if (SWIG_ConvertPtr($input, (void **) &temp_vec, $descriptor(mfem::Vector *), 0 |0) == -1) {
           if (!PyArray_Check($input) || !PyArray_ISFLOAT($input)){	 
              PyErr_SetString(PyExc_ValueError, "Expecting a list/const *double/Vector/numpy float array");
              return NULL;
	   } else {
  	      //std::cout << "Calling numpy data(float)\n";	     
              $1 = (double *) PyArray_DATA((PyArrayObject *)$input);	     
	   }
       } else {
	 //std::cout << "Calling Vector::GetData\n";
         $1 = temp_vec -> GetData();
       }
    } else {
      // const (
      //std::cout << "Using const double\n";      
      $1 = temp_ptr;
   }
  } else {
    //std::cout << "Using List\n";    
    int l = PyList_Size($input);
    $1 = (double *) malloc((l)*sizeof(double));
    for (i = 0; i < l; i++) {
      PyObject *s = PyList_GetItem($input,i);
      if (PyInt_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
      } else if (PyFloat_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
      } else {
        free($1);      
        PyErr_SetString(PyExc_ValueError, "List items must be integer/float");
        return NULL;
      }
    }
  }
}
%typemap(typecheck) (TYPENAME) {
   void *ptr;
   if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(const double *), 0 |0) == -1) { 
      PyErr_Clear();
      if (!PyList_Check($input)) {
         PyErr_Clear();
         if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(mfem::Vector *), 0 |0) == -1) {
            PyErr_Clear();
            if (!PyArray_Check($input) || !PyArray_ISFLOAT($input)){
              $1 = 0;	      
	    } else {
              $1 = 1;  // accept numpy float array
	    }
	 } else {
	   $1 = 1;  // accept vector
	 }
      } else {
	$1 = 1;  // acccept list
      }
   } else {
     $1 = 1;     // accept const double*
  }
}
%enddef
