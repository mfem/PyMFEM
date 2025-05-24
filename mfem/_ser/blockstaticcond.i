%module(package="mfem._ser") blockstaticcond
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/dpg/util/blockstaticcond.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_BLOCKSTATICCOND

%init %{
import_array();
%}

%inline %{
#include "miniapps/dpg/util/blockstaticcond.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "array.i"
%import "vector.i"
%import "densemat.i"
%import "operators.i"
%import "blockmatrix.i"
%import "blockoperator.i"
%import "fespace.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)


%include "miniapps/dpg/util/blockstaticcond.hpp"

#endif
