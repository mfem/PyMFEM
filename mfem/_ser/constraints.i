/*

   constraints.i

*/
%module(package="mfem._ser") constraints
%feature("autodoc", "1");
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pyintrules.hpp"
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
