//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") estimators

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
  %}
// initialization required to return numpy array from SWIG
%init %{
import_array1(-1);
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "fespace.i"
%import "bilinearform.i"
%import "gridfunc.i"
%import "../common/exception.i"

%ignore mfem::ZienkiewiczZhuEstimator::ZienkiewiczZhuEstimator(BilinearFormIntegrator &integ, GridFunction &sol, FiniteElementSpace *flux_fes);
%ignore mfem::ZienkiewiczZhuEstimator::ZienkiewiczZhuEstimator(BilinearFormIntegrator &integ, GridFunction &sol, FiniteElementSpace &flux_fes);

namespace mfem{
  %pythonprepend ZienkiewiczZhuEstimator::ZienkiewiczZhuEstimator %{
     if own_flux_fes: flux_fes.thisown=0
  %}
}

/* this is to ignroe final keyword */
# define final
%include "fem/estimators.hpp"


namespace mfem{
  %extend ZienkiewiczZhuEstimator{
     ZienkiewiczZhuEstimator(mfem::BilinearFormIntegrator &integ,
  			     mfem::GridFunction &sol,
  			     mfem::FiniteElementSpace *flux_fes,
			     bool own_flux_fes = false){
       if (own_flux_fes){
           return new mfem::ZienkiewiczZhuEstimator(integ, sol, flux_fes);
       } else {
           return new mfem::ZienkiewiczZhuEstimator(integ, sol, *flux_fes);
       }
     }
  };
}

