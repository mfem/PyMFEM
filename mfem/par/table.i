%module table
%{
#include "iostream_typemap.hpp"    
#include "general/table.hpp"
#include <iostream>
%}

%import array.i
%import "ostream_typemap.i"

// recipe for ostream
%typemap(in) std::ostream& (boost_ofdstream *stream=NULL) {
  FILE *f=PyFile_AsFile($input); // Verify the semantics of this
  if (!f) {
    SWIG_Error(SWIG_TypeError, "File object expected.");
    SWIG_fail;
  }
  else {
    // If threaded incrment the use count
    stream = new boost_ofdstream(fileno(f), io::never_close_handle);
    $1 = new std::ostream(stream);
  }
}
%typemap(freearg) std::ostream& {
  delete $1;
  delete stream$argnum;
}

%include "general/table.hpp"
