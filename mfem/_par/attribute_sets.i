//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") attribute_sets
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_MESH_ATTRIBUTE_SETS

%init %{
import_array1(-1);
%}

%include "exception.i"
%include "../common/exception.i"

%import "array.i"
%import "arrays_by_name.i"
%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

//
//  AttributeSets::GetAttributeSetNames returns Python set
//
%typemap(out) std::set<std::string>{
    resultobj = PySet_New(NULL);
    for (auto const &set_name : *(& $1)){
      std::ostringstream name2;
      name2 << set_name;
      PySet_Add(resultobj, PyUnicode_FromString(name2.str().c_str()));
    }
}

%include "mesh/attribute_sets.hpp"

#endif
