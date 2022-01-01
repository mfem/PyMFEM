%module(package="mfem._ser") element

%{
#include <iostream>  
#include "mfem.hpp"
#include "pyoperator.hpp"      
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"

%import "globals.i"
%import "array.i"
%import "densemat.i"
%import "geom.i"
%import "table.i"
%import "hash.i"
%import "../common/exception.i"

%include "../common/deprecation.i"
DEPRECATED_OVERLOADED_METHOD(mfem::Element::GetNFaces,
    	                     Element::GetNFaces(int & nFaceVertices) is deprecated,
			     len(args) == 1)

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
