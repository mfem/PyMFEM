//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") tetrahedron
%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}
%include "exception.i"
%import "fe.i"
%import "element.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

LIST_TO_INTARRAY_IN(const int *ind, 2)
INTARRAY_OUT_TO_TUPLE(int *GetVertices, 2)

%include "../common/deprecation.i"
DEPRECATED_OVERLOADED_METHOD(mfem::Tetrahedron::GetNFaces,
    	                     Tetrahedron::GetNFaces(int & nFaceVertices) is deprecated,
			     len(args) == 1)

%include "mesh/tetrahedron.hpp"
