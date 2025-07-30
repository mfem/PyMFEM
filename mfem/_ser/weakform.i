%module(package="mfem._ser") weakform
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/dpg/util/weakform.hpp"
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
#include "miniapps/dpg/util/weakform.cpp"
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
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)

%import "../common/object_array_typemap.i"
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementSpace*>&,
			   FiniteElementSpace*)
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementCollection*>&,
			   FiniteElementCollection*)


%pythonprepend mfem::DPGWeakForm::DPGWeakForm %{
  if len(args) > 0:
     fes_, fecol_ = args
     self._fes = fes_
     self._fecol = fecol_
%}

%pythonprepend mfem::DPGWeakForm::SetSpaces %{
  self._fes = fes_
  self._fecol = fecol_
%}


%include "miniapps/dpg/util/weakform.hpp"

#endif
