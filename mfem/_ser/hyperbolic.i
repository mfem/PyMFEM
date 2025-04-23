%module(package="mfem._ser") hyperbolic
%feature("autodoc", "1");

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pynonlininteg.hpp"  
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_HYPERBOLIC

%init %{
import_array();
%}

%include "exception.i"
%include "std_string.i"
%include "../common/exception.i"

%import "array.i"
%import "vector.i"
%import "densemat.i"
%import "eltrans.i"
%import "nonlininteg.i"

%include "fem/hyperbolic.hpp"

#endif



