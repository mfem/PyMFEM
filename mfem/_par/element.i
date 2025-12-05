//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") element

%{
#include <iostream>
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"

%import "globals.i"
%import "array.i"
%import "densemat.i"
%import "geom.i"
%import "table.i"
%import "hash.i"
%import "../common/exception.i"

%include "mesh/element.hpp"

%extend mfem::Element {
  PyObject* GetVerticesArray(void) const{
     int L = self->GetNVertices();
     int n;
     const int *v = self->GetVertices();
     npy_intp dims[] = {L};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_INT);
     int *x    = (int*)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
     for (n = 0; n < L; n++) {
        x[n] = v[n];
     }
     return array;
  }
};
