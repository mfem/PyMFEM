/*

   complex_fem.i

*/
%module(package="mfem._ser") complex_fem
%feature("autodoc", "1");
%{
#include "fem/complex_fem.hpp"  
#include "linalg/complex_operator.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "pycoefficient.hpp"  
  %}
%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "fespace.i"
%import "coefficient.i"
%import "complex_operator.i"
%import "bilinearform.i"
%import "operators.i"
%import "sparsemat.i"

%include "../common/complex_fem_ext.i"
%include "fem/complex_fem.hpp"
