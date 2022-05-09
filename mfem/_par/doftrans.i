%module(package="mfem._par") doftrans
%{
#include  "mfem.hpp"
#include "pyoperator.hpp"      
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%import "fe.i"
%import "intrules.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"


%include "fem/doftrans.hpp"
