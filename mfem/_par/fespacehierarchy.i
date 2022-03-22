/*

   fespacehierarchy.i

*/
%module(package="mfem._par") fespacehierarchy
%feature("autodoc", "1");
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}
%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "bilinearform.i"

%include "fem/fespacehierarchy.hpp"
