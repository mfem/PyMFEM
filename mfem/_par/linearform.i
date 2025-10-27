//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par", directors="1") linearform
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pylininteg.hpp"
%}

%init %{
import_array1(-1);
%}

//%include "../common/cpointers.i"
//%import "cpointers.i"
%import "coefficient.i"
%import "array.i"
%import "mesh.i"
%import "intrules.i"
%import "fe.i"
%import "fe_coll.i"
%import "densemat.i"
%import "sparsemat.i"
%import "vector.i"
%import "eltrans.i"
%import "lininteg.i"
 //%import "fem/fespace.hpp"

%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }
}

 //%include "fem/coefficient.hpp"
namespace mfem {
%pythonprepend LinearForm::AddDomainIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    lfi = args[0]
    self._integrators.append(lfi)
    lfi.thisown=0
   %}
%pythonprepend LinearForm::AddBoundaryIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    lfi = args[0]
    self._integrators.append(lfi)
    lfi.thisown=0
   %}
%pythonprepend LinearForm::AddBdrFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    lfi = args[0]
    self._integrators.append(lfi)
    lfi.thisown=0
   %}
%pythonprepend LinearForm::AddInteriorFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(lfi)
    lfi.thisown=0
   %}
}
%include "fem/linearform.hpp"


