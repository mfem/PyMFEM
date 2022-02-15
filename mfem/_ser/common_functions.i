%module(package="mfem._ser") common_functions

%{
#include <iostream>
#include <fstream>
#include "mfem.hpp"
#include "../common/io_stream.hpp"          
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"     
  %}
%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%import "globals.i"
%import "array.i"
%import "globals.i"
%ignore "";
%rename("%s") mfem;
%rename("%s") mfem::Add;
%rename("%s") mfem::Transpose;
%rename("%s") mfem::Mult;
%rename("%s") mfem::RAP;
%rename("%s") mfem::InnerProduct;
%rename("%s") mfem::OuterProduct;
%rename("%s") mfem::TransposeMult;
%rename("%s") mfem::Mult_AtDA;
%newobject mfem::Add;
%newobject mfem::RAP;
%include "general/table.hpp"
%include "linalg/vector.hpp"
%include "linalg/densemat.hpp"
%include "linalg/sparsemat.hpp"
%include "linalg/blockmatrix.hpp"
