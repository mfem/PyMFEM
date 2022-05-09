%module (package="mfem._par") eltrans

%{
#include "mfem.hpp"
#include "pyoperator.hpp"  
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"  
%}

%init %{
import_array();
%}

%include "exception.i"

%import "globals.i"
%import "array.i"
%import "vector.i"
%import "densemat.i"
%import "fe.i"
%import "intrules.i"
%import "geom.i"

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

%include "fem/eltrans.hpp"

