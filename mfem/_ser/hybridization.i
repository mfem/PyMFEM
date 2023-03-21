%module(package="mfem._ser")  hybridization
%{
#include "mfem.hpp"  
#include "pyoperator.hpp"        
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"
%include "exception.i"
%import "fespace.i"
%import "bilininteg.i"
%import "../common/exception_director.i"

%include "fem/hybridization.hpp"
