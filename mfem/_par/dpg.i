%module(package="mfem._par") dpg
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pylininteg.hpp"    
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_WEAKFORM


%init %{
import_array();
%}

%inline %{
#include "miniapps/dpg/util/pml.cpp"
#include "miniapps/dpg/util/blockstaticcond.cpp"
#include "miniapps/dpg/util/complexstaticcond.cpp"    
#include "miniapps/dpg/util/weakform.cpp"
#include "miniapps/dpg/util/complexweakform.cpp"  
#include "miniapps/dpg/util/pweakform.cpp"
#include "miniapps/dpg/util/pcomplexweakform.cpp"    
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "array.i"
%import "vector.i"
%import "blockvector.i"
%import "bilininteg.i"
%import "lininteg.i"
%import "fespace.i"
%import "blockmatrix.i"
%import "operators.i"
%import "blockmatrix.i"
%import "blockoperator.i"
%import "complex_densemat.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)

%ignore mfem::BlockStaticCondensation::ConvertListToReducedTrueDofs;
%ignore mfem::ComplexBlockStaticCondensation::ConvertListToReducedTrueDofs;

%include "miniapps/dpg/util/pml.hpp"
%include "miniapps/dpg/util/blockstaticcond.hpp"
%include "miniapps/dpg/util/complexstaticcond.hpp"
%include "miniapps/dpg/util/weakform.hpp"
%include "miniapps/dpg/util/complexweakform.hpp"
%include "miniapps/dpg/util/pweakform.hpp"
%include "miniapps/dpg/util/pcomplexweakform.hpp"


#endif
