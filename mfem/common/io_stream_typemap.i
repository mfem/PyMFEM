//usage OSTREAM_TYPEMAP(std::ostream&)
%define OSTREAM_TYPEMAP(T)
%typemap(in) T (PyMFEM::wFILE *temp  = 0, std::ofstream  out){
  if (SWIG_ConvertPtr($input, (void **) &temp, $descriptor(PyMFEM::wFILE *), 0 | 0) == -1) {
     SWIG_exception(SWIG_ValueError,"io_stream object is expected.");      
     return NULL;
  }  

  if (temp->isSTDOUT() == 1) {
      $1 = &std::cout;
  }
  else {
    out.open(temp->getFilename());
    out.precision(temp->getPrecision());
    $1 = &out;
  }
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_STRING_ARRAY) std::ostream& {
  void *ptr;
  if (SWIG_ConvertPtr($input, (void **) &ptr, $descriptor(PyMFEM::wFILE *), 0 |0) == -1) {       
    PyErr_Clear();
    $1 = 0;
  } else {
    $1 = 1;    
  }
}
%typemap(freearg) std::ostream& {
  if (temp$argnum) {
     if (temp$argnum->isSTDOUT() != 1) {
         out$argnum.close();
     }
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

