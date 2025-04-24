%module(package="mfem._par") blockstaticcond
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
