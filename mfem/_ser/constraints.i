/*

   constraints.i

*/
%module(package="mfem._ser") constraints
%feature("autodoc", "1");
%{
#include "linalg/constraints.hpp"
#include "fem/restriction.hpp"
#include "fem/linearform.hpp"      
#include "linalg/complex_operator.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pysolvers.hpp"
  
%}
%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "fespace.i"
%import "operators.i"
%import "sparsemat.i"
%import "solvers.i"

%include "linalg/constraints.hpp"
