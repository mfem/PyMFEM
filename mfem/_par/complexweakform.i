%module(package="mfem._par") complexweakform
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/dpg/util/complexweakform.hpp"
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
#ifdef FILE_EXISTS_MINIAPPS_DPG_UTIL_COMPLEXWEAKFORM

%init %{
import_array();
%}

%inline %{
#include "miniapps/dpg/util/complexweakform.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "gridfunc.i"
%import "mesh.i"
%import "solvers.i"
%import "operators.i"
%import "blockmatrix.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)

%import "../common/object_array_typemap.i"
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::FiniteElementSpace*>&,
			   FiniteElementSpace*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::FiniteElementCollection*>&,
			   FiniteElementCollection*)


%pythonprepend mfem::ComplexDPGWeakForm::ComplexDPGWeakForm %{
  if len(args) > 0:
     fes_, fecol_ = args
     self._fes = fes_
     self._fecol = fecol_
%}

%pythonprepend mfem::ComplexDPGWeakForm::SetSpaces %{
  self._fes = fes_
  self._fecol = fecol_
%}


%include "miniapps/dpg/util/complexweakform.hpp"

#endif
