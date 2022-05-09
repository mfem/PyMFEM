%module(package="mfem._par") io_stream

%feature("autodoc", "1");

%{
#include  "fstream"
#include  "iostream"
#include  "string"
#include  "../common/io_stream.hpp"    
%}

%include "../common/io_stream.hpp"

%pythoncode %{
  STDOUT = wFILE()
%}
//
