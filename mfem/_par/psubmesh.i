//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") psubmesh

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MESH_SUBMESH_SUBMESH

%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%import "submesh.i"
%import "pmesh.i"
%import "pgridfunc.i"
%import "ptransfermap.i"

%include "mesh/submesh/psubmesh.hpp"

#endif
