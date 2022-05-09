//
// const_intptr_typemap
//
// This typemap accept list/const */Array<int>/numpy float array for int32 *  (used in mesh.i)
//
// It is asumed that data is copied inside called metood.
//
// example: CONST_INTPTR_IN(const int *vi)
//
%define CONST_INTPTR_IN(TYPENAME)
%typemap(in) (TYPENAME) (mfem::Array<int> *temp_arr, int * temp_ptr, bool is_allocated=false, bool is_tuple=false)
{
  int i;
  if (!PyList_Check($input) && !PyTuple_Check($input)) {
    if (SWIG_ConvertPtr($input, (void **) &temp_ptr, $descriptor(const int *), 0 |0) == -1) {
       if (SWIG_ConvertPtr($input, (void **) &temp_arr, $descriptor(mfem::Array<int> *), 0 |0) == -1) {
         if (!PyArray_Check($input) || !PyArray_ISINTEGER(reinterpret_cast<PyArrayObject *>($input))){
            PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple/const *int/Vector/numpy int array");
            return NULL;
	 } else {
	   //std::cout << "Calling numpy data(int)\n";	     
            $1 = (int *) PyArray_DATA((PyArrayObject *)$input);
	    //std::cout << $1[0] << " " << $1[1] << " " << $1[2] << "\n";
         }	 
       } else {
	 //std::cout << "Calling Array<int>::GetData\n";
         $1 = temp_arr -> GetData();
         //std::cout << $1[0] << " " << $1[1] << " " << $1[2] << "\n";	 
       }
    } else {
      // const (
      //std::cout << "Using const int\n";      
      $1 = temp_ptr;
      //std::cout << $1[0] << " " << $1[1] << " " << $1[2] << "\n";	       
   }
  } else {
    int l = 0;
    if (PyTuple_Check($input)) {
       is_tuple = true;
       l = PyTuple_Size($input);      
    } else {
       l = PyList_Size($input);
    }
    //$1 = (int *) malloc((l)*sizeof(int));
    $1 = new int[l];
    is_allocated = true;
    for (i = 0; i < l; i++) {
      PyObject *s = (is_tuple) ? PyTuple_GetItem($input, i) : PyList_GetItem($input,i);
      if (PyInt_Check(s)) {
        $1[i] = (int)PyInt_AsLong(s);
      } else if ((PyArray_PyIntAsInt(s) != -1) || !PyErr_Occurred()) {
        $1[i] = PyArray_PyIntAsInt(s);
      } else {    
        delete[] $1;
        PyErr_SetString(PyExc_ValueError, "List/Tuple items must be integer");
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
   if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(const int *), 0 |0) == -1) { 
      PyErr_Clear();
      if (!PyList_Check($input) && !PyTuple_Check($input)) {
         PyErr_Clear();
         if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(mfem::Array<int> *), 0 |0) == -1) {
	   if (!PyArray_Check($input) || !PyArray_ISINTEGER(reinterpret_cast<PyArrayObject *>($input))){
              $1 = 0;	      
	    } else {
              $1 = 1;  // accept numpy int array
	    }
	 } else {
	   $1 = 1;  // accept array <int>
	 }
      } else {
	$1 = 1;  // acccept list/tuple
      }
   } else {
     $1 = 1;     // accept const int*
  }
}
%enddef
