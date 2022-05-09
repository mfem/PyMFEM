%module(package="mfem._par") hybridization
%{
#include "mfem.hpp"
#include "pyoperator.hpp"  
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"  
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
