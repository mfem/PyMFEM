//
//
//  ostream
//
//

//usage OSTREAM_TYPEMAP(std::ostream&)
%define OSTREAM_TYPEMAP(T)
  %typemap(in) T (PyMFEM::wFILE *temp=0, std::ofstream out_txt, mfem::ofgzstream *out_gz=0,
		  PyObject *string_io=0, std::ostringstream *stream=0,
		  PyObject *ret=0){
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
	 string_io=$input;
	 stream = new std::ostringstream();
	 int prec = 16;
	 if (PyObject_HasAttrString($input, "precision")){
	   PyObject *attr = PyObject_GetAttrString($input, "precision");
	   prec = (int)PyLong_AsLong(attr);
	   //std::cout << "setting prec" << prec << "\n";
	 }
         stream->precision(prec);
	 
      } else {
 	 // if it is string, extract filename as char*
         PyObject* str = PyUnicode_AsEncodedString($input, "utf-8", "~E~");	
         const char* filename = PyBytes_AsString(str);
         temp = new PyMFEM::wFILE(filename, 16, true);
         Py_DECREF(str);	 
      }
   }

   if (stream == 0){
      if (temp->isSTDOUT() == 1) {
         $1 = &std::cout;
      } else if (temp->isGZ()){
 	 out_gz = new mfem::ofgzstream(temp->getFilename(), true);
         $1 = out_gz;	     
      } else {
         out_txt.open(temp->getFilename());
         out_txt.precision(temp->getPrecision());
         $1 = &out_txt;
      }
   } else {
      $1 = stream;
   }
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_STRING) T {
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
        $1 = 1;
      }
    } else {
        $1 = 1;
    }
}

%typemap(freearg) T {
  if (!stream$argnum) {  
     if (temp$argnum) {
         if (temp$argnum->isSTDOUT() != 1) {
	   if (out_txt$argnum.is_open()){
	     out_txt$argnum.close();
	   }
	   if (out_gz$argnum){
	     delete out_gz$argnum;
	   }
         }
         if (temp$argnum->isTemporary()){
            delete temp$argnum;
         }
     }
  }
 }
%typemap(argout) T {
  if (stream$argnum) {  
    std::string str =  stream$argnum->str();
    const char* s = str.c_str();
    const int n = str.length();
    ret$argnum = PyObject_CallMethod(string_io$argnum, "write", "s#",
			s, static_cast<Py_ssize_t>(n));
    if (PyErr_Occurred()) {
       PyErr_SetString(PyExc_RuntimeError, "Error occured when writing IOString");
       return NULL;
    }
    delete stream$argnum;
    Py_XDECREF($result);   /* Blow away any previous result */
    $result = ret$argnum;    
  }
}
%enddef

//This macro extend class to write file and stdout (no argument)
%define OSTREAM_ADD_DEFAULT_STDOUT_FILE(class, method)
%extend mfem::class {
void method(const char *file, int precision=16){
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
void method ## GZ(const char *file, int precision=16){
  mfem::ofgzstream *ofile = new mfem::ofgzstream(file, true);
  
  if (!ofile)
     {
        std::cerr << "\nCan not produce output file: " << file << '\n' << std::endl;
        return;
      }
  ofile -> precision(precision);  
  self -> method(*ofile);
  delete ofile;
  }
void method(void){
  self -> method(std::cout);
  }
 
};
%enddef

%define OSTREAM_ADD_DEFAULT_FILE(class, method)
%extend mfem::class {
void method(const char *file, int precision=16){
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
void method ## GZ(const char *file, int precision=16){
  mfem::ofgzstream *ofile = new mfem::ofgzstream(file, true);
  if (!ofile)
     {
        std::cerr << "\nCan not produce output file: " << file << '\n' << std::endl;
        return;
      }
  ofile ->precision(precision);  
  self -> method(*ofile);
  delete ofile;
  }
};
%enddef

//
//
//  istream
//
//
//usage ISTREAM_TYPEMAP(std::istream&)
%define ISTREAM_TYPEMAP(T)
  %typemap(in) T (PyMFEM::wFILE *temp=0, std::ifstream in_txt, mfem::ifgzstream *in_gz=0,
		  std::istringstream *stream=0, Py_ssize_t len = 0){
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

	 PyObject *input_str = PyObject_CallMethod($input, "getvalue", NULL);
         if (PyErr_Occurred()) {
             PyErr_SetString(PyExc_RuntimeError, "Can not read from StringIO");
	     return NULL;
         }
	 
         char *buf = nullptr;
         PyObject *str = PyUnicode_AsUTF8String(input_str);	 
         PyBytes_AsStringAndSize(str, &buf, &len);
         stream = new std::istringstream(buf);
         Py_DECREF(str);
         Py_DECREF(input_str);	 
      } else {
 	 // if it is string, extract filename as char*
         PyObject* str = PyUnicode_AsEncodedString($input, "utf-8", "~E~");	
         const char* filename = PyBytes_AsString(str);
         temp = new PyMFEM::wFILE(filename, 8, true);
         Py_DECREF(str);
      }
   }
   if (stream == 0){
     /*
      if (temp->isGZ()){
  	 in_gz = new mfem::ifgzstream(temp->getFilename());
         $1 = in_gz;
      } else {
  	 in_txt.open(temp->getFilename(), std::ifstream::in);
         in_txt.precision(temp->getPrecision());
         $1 = &in_txt;
      }
     */
      /* this will auto-detect the input file type */
      in_gz = new mfem::ifgzstream(temp->getFilename());
      $1 = in_gz;
     
      if (temp->isTemporary()){
         delete temp;
      }
   } else {
      $1 = stream;
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
        $1 = 1;
      }
    } else {
        $1 = 1;
    }
}

%typemap(freearg) T {
  if (!stream$argnum) {    
     if (temp$argnum) {
         in_txt$argnum.close();
     }
     if (in_gz$argnum){
        delete in_gz$argnum;
     }
  }  
 }
/*
%typemap(argout) T {
  if (stream$argnum) {  
    ret$argnum = PyLong_FromSsize_t(len$argnum);
    if (PyErr_Occurred()) {
       PyErr_SetString(PyExc_RuntimeError, "Error occured when writing IOString");
       return NULL;
    }
    delete stream$argnum;    
    Py_XDECREF($result);   // Blow away any previous result
    $result = ret$argnum;    
  }
}
*/
%enddef

