%module eltrans

%{
#include "general/array.hpp"
#include "fem/intrules.hpp"
#include "fem/eltrans.hpp"
#include "numpy/arrayobject.h"      
%}

%init %{
import_array();
%}

%import array.i
%import vector.i
%import densemat.i
%import fe.i
%import intrules.i

%feature("shadow") mfem::ElementTransformation::Transform %{
def Transform(self, *args):
    from .vector import Vector
    from .intrules import IntegrationPoint
    if isinstance(args[0], IntegrationPoint):
        vec = Vector()
        _eltrans.ElementTransformation_Transform(self, args[0], vec)
        vec.thisown = 0      
        return vec.GetDataArray()
    else:
        return _eltrans.ElementTransformation_Transform(self, *args)
%}

%include "fem/eltrans.hpp"

