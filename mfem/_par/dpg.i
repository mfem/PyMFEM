%module(package="mfem._par") dpg

%feature("autodoc", "1");

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
using namespace mfem;
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_WEAKFORM


%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%inline %{
#include "miniapps/dpg/util/pml.cpp"
#include "miniapps/dpg/util/blockstaticcond.cpp"
#include "miniapps/dpg/util/complexstaticcond.cpp"
#include "miniapps/dpg/util/weakform.cpp"
#include "miniapps/dpg/util/complexweakform.cpp"
#include "miniapps/dpg/util/pweakform.cpp"
#include "miniapps/dpg/util/pcomplexweakform.cpp"
%}

%include "../common/dpg_common.i"

%include "miniapps/dpg/util/pml.hpp"
%include "miniapps/dpg/util/blockstaticcond.hpp"
%include "miniapps/dpg/util/complexstaticcond.hpp"
%include "miniapps/dpg/util/weakform.hpp"
%include "miniapps/dpg/util/complexweakform.hpp"
%include "miniapps/dpg/util/pweakform.hpp"
%include "miniapps/dpg/util/pcomplexweakform.hpp"


#endif
