%module(package="mfem._par") istream_typemap

//recipe for istream
%typemap(in) std::istream& (boost_ifdstream *stream=NULL) {
  FILE *f=PyFile_AsFile($input); // Verify that this returns NULL for non-files
  if (!f) {
    SWIG_Error(SWIG_TypeError, "File object expected.");  
    SWIG_fail;
  }
  else {
    stream = new boost_ifdstream(fileno(f), io::never_close_handle);
    $1 = new std::istream(stream);
  }
}
%typemap(typecheck) std::istream& {
  if (PyFile_Check($input)){
    $1 = 1;
  } else {
    $1 = 0;
  }
}
%typemap(freearg) std::istream& {
  delete $1;
  delete stream$argnum;
}


