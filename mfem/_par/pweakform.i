%module(package="mfem._par") pweakform
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/dpg/util/pweakform.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_PWEAKFORM

%init %{
import_array();
%}

%inline %{
#include "miniapps/dpg/util/pweakform.cpp"
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

OSTREAM_TYPEMAP(std::ostream&)

%import "../common/object_array_typemap.i"
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::ParFiniteElementSpace*>&,
			   ParFiniteElementSpace*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::FiniteElementCollection*>&,
			   FiniteElementCollection*)


%pythonprepend mfem::ParDPGWeakForm::ParDPGWeakForm %{
    if len(args) > 0:
     trial_fes_, fecol_ = args
     self._fes = trial_pfes_
     self._fecol = fecol_

%}

%pythonprepend mfem::ParDPGWeakForm::SetParSpaces %{
  self._fes = trial_pfes_
  self._fecol = fecol_
%}

%include "miniapps/dpg/util/pweakform.hpp"

#endif
