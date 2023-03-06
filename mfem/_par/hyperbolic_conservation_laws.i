%module(package="mfem._ser", directors="0")  hyperbolic_conservation_laws
%begin %{
#ifndef PY_SSIZE_T_CLEAN  
#define PY_SSIZE_T_CLEAN
#endif
%}

%{

  #include <fstream>  
  #include <iostream>
  #include <sstream>
  #include <limits>
  #include <cmath>
  #include <cstring>
  #include <ctime>
  #include "mfem/mfem.hpp"  
  #include "../common/pycoefficient.hpp"
  #include "pyoperator.hpp"  
  #include "numpy/arrayobject.h"
  #include "../common/io_stream.hpp"
  #include "fem/hyperbolic_conservation_laws.hpp"  
//  using namespace mfem;  
%}


// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%include "exception.i"
%include "std_string.i"

%import "array.i"
%import "vector.i"
%import "coefficient.i"
%import "mesh.i"
%import "fespace.i"
%import "bilininteg.i"
%import "linearform.i"
%import "fespace.i"
%import "fe_coll.i"
%import "intrules.i"
%import "densemat.i"
%import "sparsemat.i"
%import "lininteg.i"
%import "eltrans.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%include "fem/hyperbolic_conservation_laws.hpp"
