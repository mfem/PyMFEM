%module  hybridization
%{
#include "fem/linearform.hpp"
#include "fem/hybridization.hpp"
#include "numpy/arrayobject.h"      
%}

%init %{
import_array();
%}

%import "fespace.i"
%import "bilininteg.i"

%include "fem/hybridization.hpp"
