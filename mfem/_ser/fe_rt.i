%module(package="mfem._ser") fe_rt
%{
#include  "mfem.hpp"
#include "fem/fe/fe_rt.hpp"
#include "mesh/hexahedron.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%import "fe_base.i"
%import "element.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"


%include "fem/fe/fe_rt.hpp"