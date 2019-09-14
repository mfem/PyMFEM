%module(package="mfem._ser") io_stream

%feature("autodoc", "1");

%{
#include  "fstream"
#include  "iostream"
#include  "string"
#include  "io_stream.hpp"    
%}

%include "io_stream.hpp"

%pythoncode %{
  STDOUT = wFILE('__stdout__', 8)
%}
//
