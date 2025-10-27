//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module (package="mfem._ser") eltrans

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/io_stream.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"

%import "globals.i"
%import "array.i"
%import "vector.i"
%import "densemat.i"
%import "fe.i"
%import "intrules.i"
%import "geom.i"
%import "../common/exception.i"
%import "../common/mfem_config.i"

%feature("shadow") mfem::ElementTransformation::Transform %{
def Transform(self, *args):
    from .vector import Vector
    from .intrules import IntegrationPoint
    if isinstance(args[0], IntegrationPoint):
        vec = Vector()
        _eltrans.ElementTransformation_Transform(self, args[0], vec)
        ret = vec.GetDataArray().copy()
        return ret
    else:
        return _eltrans.ElementTransformation_Transform(self, *args)
%}

%include "../common/deprecation.i"
DEPRECATED_METHOD(mfem::IsoparametricTransformation::FinalizeTransformation())


%ignore mfem::ElementTransformation::TransformBack;
%ignore mfem::IsoparametricTransformation::TransformBack;

%include "../common/kernel_dispatch.i"
%include "fem/eltrans.hpp"

//
//  special handling for TransformBack (this is because tol_0 is protected)
//
namespace mfem{
  #ifdef MFEM_USE_DOUBLE
  %extend IsoparametricTransformation{
     virtual int _TransformBack(const Vector &pt, IntegrationPoint &ip,
                                const real_t phys_tol = 1e-15){
       return self-> TransformBack(pt, ip, phys_tol);
     }
   };  //end of extend
  #elif defined(MFEM_USE_SINGLE)
  %extend IsoparametricTransformation{
     virtual int _TransformBack(const Vector &pt, IntegrationPoint &ip,
                                const real_t phys_tol = 1e-7){
       return self-> TransformBack(pt, ip, phys_tol);
     }
   };  //end of extend
  #endif
 } //end of namespace

%pythoncode %{
if hasattr(IsoparametricTransformation, "_TransformBack"):
    IsoparametricTransformation.TransformBack = IsoparametricTransformation._TransformBack
%}


