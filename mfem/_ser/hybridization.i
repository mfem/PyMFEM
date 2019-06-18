%module(package="mfem._ser")  hybridization
%{
#include "fem/gridfunc.hpp"  
#include "fem/linearform.hpp"
#include "fem/hybridization.hpp"
#include "numpy/arrayobject.h"      
%}

%init %{
import_array();
%}

%include "exception.i"
%import "fespace.i"
%import "bilininteg.i"
%import "../common/exception_director.i"

%include "fem/hybridization.hpp"
