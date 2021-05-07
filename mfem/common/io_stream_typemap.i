//usage OSTREAM_TYPEMAP(std::ostream&)
%define OSTREAM_TYPEMAP(T)
  %typemap(in) T (PyMFEM::wFILE *temp=0, std::ofstream out, int use_stringio=0, PyObject*string_io=0, std::ostringstream stream, PyObject* ret=0){
   //  PyMFEM::wFILE or string argument or StringIO
   if (SWIG_ConvertPtr($input, (void **) &temp, $descriptor(PyMFEM::wFILE *), 0 | 0) == -1) {
      if (!PyString_Check($input) && !PyUnicode_Check($input)) {	
 	 // not string, check if it is StringIO
         PyObject* module = PyImport_ImportModule("io");
         if (!module){
             PyErr_SetString(PyExc_RuntimeError, "Can not load io module");
             return NULL;
         }      
         PyObject* cls = PyObject_GetAttrString(module, "StringIO");
         if (!cls){
             PyErr_SetString(PyExc_RuntimeError, "Can not load StringIO");
             return NULL;
         }      
         int check = PyObject_IsInstance($input, cls);
         Py_DECREF(module);
         if (! check){
            SWIG_exception(SWIG_ValueError,"First argument must be string/wFILE/IOString");
            return NULL;
         }
         use_stringio=1;
	 string_io=$input;
      } else {
 	 // if it is string, extract filename as char*
         PyObject* str = PyUnicode_AsEncodedString($input, "utf-8", "~E~");	
         const char* filename = PyBytes_AsString(str);
         temp = new PyMFEM::wFILE(filename, 8, true);	
      }
   }

   if (use_stringio == 0){
      if (temp->isSTDOUT() == 1) {
         $1 = &std::cout;
      } else {
         out.open(temp->getFilename());
         out.precision(temp->getPrecision());
         if (temp->isTemporary()){
            delete temp;
         }
         $1 = &out;
      }
   } else {
      $1 = &stream;
   }
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_STRING_ARRAY) T {
  void *ptr;
  //std::string *ptr2 = (std::string *)0;
  if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(PyMFEM::wFILE *), 0 |0) == -1) {
      PyErr_Clear();
      if (!PyString_Check($input) && !PyUnicode_Check($input)) {	
 	 // not string
         $1 = 1;	   	
         PyObject* module = PyImport_ImportModule("io");
         if (!module){
             $1 = 0;	   
         }      
         PyObject* cls = PyObject_GetAttrString(module, "StringIO");
         if (!cls){
             $1 = 0;	   	   
         }      
         int check = PyObject_IsInstance($input, cls);
         Py_DECREF(module);
         if (! check){
             $1 = 0;	   	   	   
         }
      } else {
        $1 = 0;
      }
    } else {
        $1 = 1;
    }
}

%typemap(freearg) T {
  if (use_stringio$argnum == 0) {  
     if (temp$argnum) {
         if (temp$argnum->isSTDOUT() != 1) {
            out$argnum.close();
         }
     }
  }
 }

%typemap(argout) T {
  if (use_stringio$argnum == 1) {  
    std::string str =  stream$argnum.str();
    const char* s = str.c_str();
    const int n = str.length();
    ret$argnum = PyObject_CallMethod(string_io$argnum, "write", "s#",
			s, static_cast<Py_ssize_t>(n));
    if (PyErr_Occurred()) {
       PyErr_SetString(PyExc_RuntimeError, "Error occured when writing IOString");
       return NULL;
    }
    Py_XDECREF($result);   /* Blow away any previous result */
    $result = ret$argnum;    
  }
}
%enddef

//This macro extend class to write file and stdout (no argument)
%define OSTREAM_ADD_DEFAULT_STDOUT_FILE(class, method)
%extend mfem::class {
void method(const char *file, int precision=8){
  std::ofstream ofile(file);
  if (!ofile)
     {
        std::cerr << "\nCan not produce output file: " << file << '\n' << std::endl;
        return;
      }
  ofile.precision(precision);    
  self -> method(ofile);
  ofile.close();
  }
void method(void){
  self -> method(std::cout);
  }
 
};
%enddef

%define OSTREAM_ADD_DEFAULT_FILE(class, method)
%extend mfem::class {
void method(const char *file, int precision=8){
  std::ofstream ofile(file);
  if (!ofile)
     {
        std::cerr << "\nCan not produce output file: " << file << '\n' << std::endl;
        return;
      }
  ofile.precision(precision);  
  self -> method(ofile);
  ofile.close();
  }
};
%enddef


