%module(package="mfem._ser") fe_coll
%{
#include "fem/fem.hpp"
#include "fem/fe_coll.hpp"
#include "fem/fespace.hpp"
#include "fem/eltrans.hpp"
#include "fem/coefficient.hpp"
#include "fem/intrules.hpp"  
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "numpy/arrayobject.h"      
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

 //%inline %{
 //  typedef mfem::L2_FECollection mfem::DG_FECollection;
 // %}

%include "fem/fe_coll.hpp"
%pythoncode %{
  DG_FECollection = L2_FECollection
%}   

