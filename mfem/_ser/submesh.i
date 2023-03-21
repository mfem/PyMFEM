%module(package="mfem._ser") submesh

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp" 
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MESH_SUBMESH_SUBMESH

%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%import "mesh.i"
%import "transfermap.i"

%include "mesh/submesh/submesh.hpp"

#endif

