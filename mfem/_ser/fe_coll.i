%module(package="mfem._ser") fe_coll
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}

%include "exception.i"
%import "mesh.i"
%import "array.i"
%import "matrix.i"
%import "intrules.i"
%import "coefficient.i"
%import "fe.i"
%import "densemat.i"
%import "sparsemat.i"
%import "vector.i"
%import "eltrans.i"
%import "lininteg.i"
%import "../common/exception.i"

/* define FiniteElementCollectionPtrArray */
%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT_SWIGOBJ(mfem::FiniteElementCollection *, 1)
%import "../common/array_instantiation_macro.i"
IGNORE_ARRAY_METHODS(mfem::FiniteElementCollection *)
INSTANTIATE_ARRAY0(FiniteElementCollection *, FiniteElementCollection, 1)

%include "fem/fe_coll.hpp"
%pythoncode %{
  DG_FECollection = L2_FECollection
%}

