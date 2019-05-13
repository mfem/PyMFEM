%module(package="mfem._par") hybridization
%{
#include "linalg/hypre.hpp"
#include "linalg/handle.hpp"  
#include "fem/gridfunc.hpp"  
#include "fem/linearform.hpp"
#include "fem/hybridization.hpp"
#include "numpy/arrayobject.h"
  
%}

%init %{
import_array();
%}

%import "handle.i"
%import "fespace.i"
%import "bilininteg.i"
%import "hypre.i"

%pointer_class(int, intp);

%include "fem/hybridization.hpp"
