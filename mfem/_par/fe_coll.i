%module(package="mfem._par") fe_coll
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "mfem.hpp"  
#include "pyoperator.hpp"  
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}

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

%pointer_class(int, intp);

%include "fem/fe_coll.hpp"
%pythoncode %{
  DG_FECollection = L2_FECollection
%}   
