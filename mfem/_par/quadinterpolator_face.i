//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") quadinterpolator_face
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_FEM_QUADINTERPOLATOR_FACE

%init %{
import_array1(-1);
%}
%include "exception.i"
%import "fe.i"
%import "fe_fixed_order.i"
%import "element.i"
%import "mesh.i"
%import "qspace.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "../common/numpy_int_typemap.i"

%include "fem/quadinterpolator_face.hpp"

#endif
