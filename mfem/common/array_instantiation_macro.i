%define INSTANTIATE_ARRAY0(XXX, YYY)
%template(##YYY##Array) mfem::Array<mfem::XXX>;
%extend mfem::Array<mfem::XXX> { 
  PyObject * __getitem__(PyObject* param) {
    int len = self->Size();    
    if (PySlice_Check(param)) {
        long start = 0, stop = 0, step = 0, slicelength = 0;
        int check;

	//%#ifdef TARGET_PY3
   	check = PySlice_GetIndicesEx(param, len, &start, &stop, &step,
				     &slicelength);
        //%#else
   	//check = PySlice_GetIndicesEx((PySliceObject*)param, len, &start, &stop, &step,
	//			     &slicelength);
	//%#endif

	if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<T> failed.");
            return NULL; 
	}
	if (step == 1) {
	  mfem::Array<mfem::XXX> *vec;
	  vec = new mfem::Array<mfem::XXX>(self->GetData() +  start, slicelength);
	  return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<mfem::XXX> *), 1);  
	} else {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<T> with stride>1 not supported.");
	    return NULL;
	}
    } else {
        PyErr_Clear();
        long idx = PyInt_AsLong(param);
        if (PyErr_Occurred()) {
           PyErr_SetString(PyExc_ValueError, "Argument must be either int or slice");
            return NULL; 	
        }
        int own =  (self -> OwnsData()) ? 0 : 1;
        if (idx >= 0){
          return SWIG_NewPointerObj(SWIG_as_voidptr(&(self->operator[](idx))), $descriptor(XXX *), own);
        } else {
	  return SWIG_NewPointerObj(SWIG_as_voidptr(&(self->operator[](len+idx))), $descriptor(XXX *), own);
	}
    }
  }
 };
%enddef

%define INSTANTIATE_ARRAY(XXX)
INSTANTIATE_ARRAY0(XXX, XXX)  
%enddef

%define INSTANTIATE_ARRAY_INT
%template(intArray) mfem::Array<int>;
%extend mfem::Array<int> { 
  PyObject * __getitem__(PyObject* param) {
    int len = self->Size();    
    if (PySlice_Check(param)) {
        long start = 0, stop = 0, step = 0, slicelength = 0;
        int check;

	//%#ifdef TARGET_PY3
   	check = PySlice_GetIndicesEx(param, len, &start, &stop, &step,
				     &slicelength);
        //%#else
   	//check = PySlice_GetIndicesEx((PySliceObject*)param, len, &start, &stop, &step,
        //				     &slicelength);
        //%#endif

	if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<T> failed.");
            return NULL; 
	}
	if (step == 1) {
            mfem::Array<int> *vec;
            vec = new mfem::Array<int>(self->GetData() +  start, slicelength);
            return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<int> *), 1);  
	} else {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<T> with stride>1 not supported.");
	    return NULL;
	}
    } else {
        PyErr_Clear();
        long idx = PyInt_AsLong(param);
        if (PyErr_Occurred()) {
           PyErr_SetString(PyExc_ValueError, "Argument must be either int or slice");
            return NULL; 	
        }
        if (idx >= 0){
          return PyLong_FromLong(self->operator[](idx));
        } else {
          return PyLong_FromLong(self->operator[](len+idx));
	}
    }
  }
 };
%enddef

%define INSTANTIATE_ARRAY_DOUBLE
%template(doubleArray) mfem::Array<double>;
%extend mfem::Array<double> { 
  PyObject * __getitem__(PyObject* param) {
    int len = self->Size();    
    if (PySlice_Check(param)) {
        long start = 0, stop = 0, step = 0, slicelength = 0;
        int check;

	//%#ifdef TARGET_PY3
   	check = PySlice_GetIndicesEx(param, len, &start, &stop, &step,
				     &slicelength);
        //%#else
   	//check = PySlice_GetIndicesEx((PySliceObject*)param, len, &start, &stop, &step,
	//			     &slicelength);
	//%#endif

	if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<T> failed.");
            return NULL; 
	}
	if (step == 1) {
            mfem::Array<double> *vec;
            vec = new mfem::Array<double>(self->GetData() +  start, slicelength);
            return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<double> *), 1);  
	} else {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<T> with stride>1 not supported.");
	    return NULL;
	}
    } else {
        PyErr_Clear();
        long idx = PyInt_AsLong(param);
        if (PyErr_Occurred()) {
           PyErr_SetString(PyExc_ValueError, "Argument must be either int or slice");
            return NULL; 	
        }
        if (idx >= 0){
          return PyFloat_FromDouble(self->operator[](idx));
        } else {
          return PyFloat_FromDouble(self->operator[](len+idx));
	}
    }
  }
 };
%enddef
