%module(package="mfem._ser") pml
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/dpg/util/pcomplexweakform.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_PML

%init %{
import_array();
%}

%inline %{
#include "miniapps/dpg/util/pml.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "array.i"
%import "vector.i"
%import "gridfunc.i"
%import "mesh.i"
%import "solvers.i"
%import "operators.i"
%import "blockmatrix.i"
%import "complex_densemat.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

%include "miniapps/dpg/util/pml.hpp"

#endif
