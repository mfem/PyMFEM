/*

   fespacehierarchy.i

*/
%module(package="mfem._par") fespacehierarchy
%feature("autodoc", "1");
%{
#include "linalg/operator.hpp"
#include "linalg/handle.hpp"
#include "fem/linearform.hpp"
#include "fem/bilinearform.hpp"    
#include "fem/fespacehierarchy.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"  
%}
%init %{
import_array();
%}

%include "exception.i"
%import "vector.i"
%import "bilinearform.i"

%include "fem/fespacehierarchy.hpp"
