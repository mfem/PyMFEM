%define INSTANTIATE_ARRAY2(XXX, YYY, ZZZ, USEPTR)
#if USEPTR == 1
%template(##ZZZ##Ptr##Array) mfem::Array<mfem::XXX>;
#else
%template(##ZZZ##Array) mfem::Array<mfem::XXX>;
#endif
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
        //                           &slicelength);
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
        int own =  0;
        swig_type_info *ty = $descriptor(const mfem::YYY *);
        #if USEPTR == 1
        if (idx >= 0){
          return SWIG_NewPointerObj(SWIG_as_voidptr((self->operator[](idx))), ty, own);
        } else {
          return SWIG_NewPointerObj(SWIG_as_voidptr((self->operator[](len+idx))), ty, own);
        }
        #else
        if (idx >= 0){
          return SWIG_NewPointerObj(SWIG_as_voidptr(&(self->operator[](idx))), ty, own);
        } else {
          return SWIG_NewPointerObj(SWIG_as_voidptr(&(self->operator[](len+idx))), ty, own);
        }
        #endif
    }
  }
 };
%enddef

%define INSTANTIATE_ARRAY(XXX)
INSTANTIATE_ARRAY0(XXX, XXX, 0)
%enddef

%define INSTANTIATE_ARRAY0(XXX, YYY, USEPTR)
INSTANTIATE_ARRAY2(XXX, YYY, YYY, USEPTR)
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
        //                                   &slicelength);
        //%#endif

        if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<int> failed.");
            return NULL;
        }
        if (step == 1) {
            mfem::Array<int> *vec;
            vec = new mfem::Array<int>(self->GetData() +  start, slicelength);
            return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<int> *), 1);
        } else {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<int> with stride>1 not supported.");
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
  PyObject* GetDataArray(void) const{
     const int * A = self->GetData();
     int L = self->Size();
     npy_intp dims[] = {L};
     if (sizeof(int) == 4){
        return  PyArray_SimpleNewFromData(1, dims, NPY_INT32, (void *)A);
     } else {
       // sizeof(int) == 8:
       return  PyArray_SimpleNewFromData(1, dims, NPY_INT64, (void *)A);
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
        //                           &slicelength);
        //%#endif

        if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<double> failed.");
            return NULL;
        }
        if (step == 1) {
            mfem::Array<double> *vec;
            vec = new mfem::Array<double>(self->GetData() +  start, slicelength);
            return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<double> *), 1);
        } else {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<double> with stride>1 not supported.");
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
  PyObject* GetDataArray(void) const{
     const double * A = self->GetData();
     int L = self->Size();
     npy_intp dims[] = {L};
     return  PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void *)A);
  }
 };
%enddef

%define INSTANTIATE_ARRAY_BOOL
%template(boolArray) mfem::Array<bool>;
%extend mfem::Array<bool> {
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
        //                           &slicelength);
        //%#endif

        if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<bool> failed.");
            return NULL;
        }
        if (step == 1) {
            mfem::Array<bool> *vec;
            vec = new mfem::Array<bool>(self->GetData() +  start, slicelength);
            return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<bool> *), 1);
        } else {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<bool> with stride>1 not supported.");
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
          if (self->operator[](idx)){
            return Py_True;
          } else {
            return Py_False;
          }
        } else {
          if (self->operator[](len+idx)){
            return Py_True;
          } else {
            return Py_False;
          }
        }
    }
  }
  PyObject* GetDataArray(void) const{
     const bool* A = self->GetData();
     int L = self->Size();
     npy_intp dims[] = {L};
     return  PyArray_SimpleNewFromData(1, dims, NPY_BOOL, (void *)A);
  }
 };
%enddef

// tool to deffine mfem::Array returining elements using numpy Scalar
 // mfem::Array(unsigned int) -> uintArray
 // NP_TYPE is things like NPY_FLOAT32
%define INSTANTIATE_ARRAY_NUMPYARRAY(XXX, YYY, NP_TYPE)
%template(##XXX##Array) mfem::Array<YYY>;
%extend mfem::Array<YYY> {
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
        //                           &slicelength);
        //%#endif

        if (check == -1) {
            PyErr_SetString(PyExc_ValueError, "Slicing mfem::Array<bool> failed.");
            return NULL;
        }
        if (step == 1) {
            mfem::Array<YYY> *vec;
            vec = new mfem::Array<YYY>(self->GetData() +  start, slicelength);
            return SWIG_NewPointerObj(SWIG_as_voidptr(vec), $descriptor(mfem::Array<YYY> *), 1);
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
        PyObject *np_val = NULL;
        PyArray_Descr *descr = NULL;
        if(! (descr = PyArray_DescrFromType(NP_TYPE))) {
           PyErr_SetString(PyExc_TypeError, "Improper descriptor");
           return NULL;
        }

        YYY *data_ptr;
        if (idx >= 0){
          data_ptr = &(self->operator[](idx));
        } else {
          data_ptr = &(self->operator[](idx+len));
        }
        np_val = PyArray_Scalar(data_ptr, descr, NULL);
       return np_val;
    }
  }
  PyObject* GetDataArray(void) const{
     const YYY * A = self->GetData();
     int L = self->Size();
     npy_intp dims[] = {L};
     return  PyArray_SimpleNewFromData(1, dims, NP_TYPE, (void *)A);
  }

 };
%enddef

%define IGNORE_ARRAY_METHODS(XXX)
%ignore mfem::Array<XXX>::Union;
%ignore mfem::Array<XXX>::Find;
%ignore mfem::Array<XXX>::FindSorted;
%ignore mfem::Array<XXX>::Sort;
%ignore mfem::Array<XXX>::DeleteFirst;
%ignore mfem::Array<XXX>::Unique;
%ignore mfem::Array<XXX>::PartialSum;
%ignore mfem::Array<XXX>::Sum;
%ignore mfem::Array<XXX>::IsSorted;
%ignore mfem::Array<XXX>::Save;
%ignore mfem::Array<XXX>::Max;
%ignore mfem::Array<XXX>::Min;
%ignore mfem::Array<XXX>::Print;
%ignore mfem::Array<XXX>::PrintGZ;
%ignore mfem::Array<XXX>::SaveGZ;
%ignore mfem::Array<XXX>::Load;
%ignore mfem::Array2D<XXX>::Print;
%ignore mfem::Array2D<XXX>::PrintGZ;
%ignore mfem::Array2D<XXX>::SaveGZ;
%enddef


%define IGNORE_ARRAY_METHODS_PREMITIVE(XXX)
%ignore mfem::Array<XXX>::PartialSum;
%ignore mfem::Array<XXX>::Sum;
%ignore mfem::Array<XXX>::IsSorted;
%ignore mfem::Array<XXX>::Save;
%ignore mfem::Array<XXX>::Max;
%ignore mfem::Array<XXX>::Min;
%ignore mfem::Array<XXX>::Print;
%ignore mfem::Array<XXX>::PrintGZ;
%ignore mfem::Array<XXX>::SaveGZ;
%ignore mfem::Array<XXX>::Load;
%ignore mfem::Array2D<XXX>::Print;
%ignore mfem::Array2D<XXX>::PrintGZ;
%ignore mfem::Array2D<XXX>::SaveGZ;
%enddef

