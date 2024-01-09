/*

   fespacehierarchy.i

*/
%module(package="mfem._ser") fespacehierarchy
%feature("autodoc", "1");
%{
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
%import "vector.i"
%import "bilinearform.i"

%include "fem/fespacehierarchy.hpp"
