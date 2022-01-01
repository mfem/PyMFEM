//
// const_doubleptr_typemap
//
// this typemap accept list/const */vector/numpy float array for double *  (used in mesh.i)
//
// It is asumed that data is copied inside called metood.
//
// example: CONST_DOUBLEPTR_IN(const double *)
//
%define CONST_DOUBLEPTR_IN(TYPENAME)
  %typemap(in) (TYPENAME) (mfem::Vector *temp_vec, double * temp_ptr, bool is_allocated=false, bool is_tuple=false)
{
  int i;
  if (!PyList_Check($input) && !PyTuple_Check($input)) {
    if (SWIG_ConvertPtr($input, (void **) &temp_ptr, $descriptor(const double *), 0 |0) == -1) {
       if (SWIG_ConvertPtr($input, (void **) &temp_vec, $descriptor(mfem::Vector *), 0 |0) == -1) {
	 if (!PyArray_Check($input) || !PyArray_ISFLOAT(reinterpret_cast<PyArrayObject *>($input))){
              PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple/const *double/Vector/numpy float array");
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
    int l = 0;
    if (PyTuple_Check($input)) {
      is_tuple = true;
      l = PyTuple_Size($input);      
    } else {
      l = PyList_Size($input);
    }      
    //std::cout << "Using List\n";    
    $1 = new double [l];
    is_allocated = true;
    for (i = 0; i < l; i++) {
      PyObject *s = (is_tuple) ? PyTuple_GetItem($input, i) : PyList_GetItem($input,i);      
      if (PyInt_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
      } else if (PyFloat_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
      } else {
        delete[] $1;      
        PyErr_SetString(PyExc_ValueError, "List/Tuple items must be integer/float");
        return NULL;
      }
    }
  }
}
%typemap(freearg) (TYPENAME)
{
 if (is_allocated$argnum)
   {
    delete[] $1;
   }
}
%typemap(typecheck) (TYPENAME) {
   void *ptr;
   if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(const double *), 0 |0) == -1) { 
      PyErr_Clear();
      if (!PyTuple_Check($input) && !PyList_Check($input)){
         PyErr_Clear();
         if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(mfem::Vector *), 0 |0) == -1) {
            PyErr_Clear();
            if (!PyArray_Check($input) || !PyArray_ISFLOAT(reinterpret_cast<PyArrayObject *>($input))){
              $1 = 0;	      
	    } else {
              $1 = 1;  // accept numpy float array
	    }
	 } else {
	   $1 = 1;  // accept vector
	 }
      } else {
	$1 = 1;  // acccept list/tuple
      }
   } else {
     $1 = 1;     // accept const double*
  }
}
%enddef
