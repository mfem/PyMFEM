%module(package="mfem._ser") gslib
%{
#include "mfem.hpp"
#include "fem/gslib.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"  
#include "../common/pycoefficient.hpp"
%}

%include "../common/mfem_config.i"

%init %{
import_array();
%}

%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"
%import vector.i
%import mesh.i
%import gridfunc.i

%include "fem/gslib.hpp"
