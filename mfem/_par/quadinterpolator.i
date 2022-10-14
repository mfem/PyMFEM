%module(package="mfem._par") quadinterpolator
%{
#include "mfem.hpp"
#include "../common/io_stream.hpp"  
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
#include "numpy/arrayobject.h"    
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_QUADINTERPOLATOR_FACE

%init %{
import_array();
%}
%include "exception.i"
%import "fe.i"
%import "fe_fixed_order.i"
%import "element.i"
%import "mesh.i"
%import "qspace.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "../common/numpy_int_typemap.i"

%include "fem/quadinterpolator.hpp"

#endif
