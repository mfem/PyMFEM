//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") complexweakform
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
import_array1(-1);
%}

%inline %{
#include "miniapps/dpg/util/complexweakform.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "fe_coll.i"
%import "gridfunc.i"
%import "mesh.i"
%import "solvers.i"
%import "operators.i"
%import "fespace.i"
%import "blockmatrix.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)

%include "../common/typemap_macros.i"
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementSpace*>&,
			 mfem::FiniteElementSpace*)
LIST_TO_MFEMOBJ_ARRAY_IN(mfem::Array<mfem::FiniteElementCollection*>&,
			 mfem::FiniteElementCollection*)


%pythonprepend mfem::ComplexDPGWeakForm::ComplexDPGWeakForm %{
  if len(args) > 0:
     fes_, fecol_ = args
     self._fes = fes_
     self._fecol = fecol_
  self._integrators = []
%}

%pythonprepend mfem::ComplexDPGWeakForm::SetSpaces %{
  self._fes = fes_
  self._fecol = fecol_
%}

%pythonappend mfem::ComplexDPGWeakForm::AddDomainLFIntegrator %{
  if lfi_r is not None:
      lfi_r.this.disown()
      self._integrators.append(lfi_r)
  if lfi_i is not None:
      lfi_i.this.disown()
      self._integrators.append(lfi_i)
%}
%pythonappend mfem::ComplexDPGWeakForm::AddTrialIntegrator %{
  if bfi_r is not None:
      bfi_r.this.disown()
      self._integrators.append(bfi_r)
  if bfi_i is not None:
      bfi_i.this.disown()
      self._integrators.append(bfi_i)
%}
%pythonappend mfem::ComplexDPGWeakForm::AddTestIntegrator %{
  if bfi_r is not None:
      bfi_r.this.disown()
      self._integrators.append(bfi_r)
  if bfi_i is not None:
      bfi_i.this.disown()
      self._integrators.append(bfi_i)
%}

%include "miniapps/dpg/util/complexweakform.hpp"

#endif
