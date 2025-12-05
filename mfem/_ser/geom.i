//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") geom
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
%}

%init %{
import_array1(-1);
%}

%include "exception.i"
%import "intrules.i"
%import "densemat.i"
%import "../common/exception.i"

%immutable RefPts;
%immutable GlobGeometryRefiner;

%import "../common/array_listtuple_typemap.i"
ARRAY_LISTTUPLE_INPUT(mfem::Geometry::Type, PyLong_AsLong)

%include "fem/geom.hpp"

%import "../common/array_listtuple_typemap.i"

namespace mfem{
   %ignore Array<Geometry::Type>::Union;
   %ignore Array<Geometry::Type>::Find;
   %ignore Array<Geometry::Type>::FindSorted;
   %ignore Array<Geometry::Type>::Sort;
   %ignore Array<Geometry::Type>::DeleteFirst;
   %ignore Array<Geometry::Type>::Unique;
   %ignore Array<Geometry::Type>::PartialSum;
   %ignore Array<Geometry::Type>::Abs;
   %ignore Array<Geometry::Type>::Sum;
   %ignore Array<Geometry::Type>::IsSorted;
   %ignore Array<Geometry::Type>::Save;
   %ignore Array<Geometry::Type>::Max;
   %ignore Array<Geometry::Type>::Min;
   %ignore Array<Geometry::Type>::IsConstant;
   %ignore Array<Geometry::Type>::Print;
   %ignore Array<Geometry::Type>::Load;
}
%template(GeometryTypeArray) mfem::Array<mfem::Geometry::Type>;
%extend mfem::Array<mfem::Geometry::Type> {
  const mfem::Geometry::Type  & __getitem__(const int i) const{
     return (* self)[i];
 }
};

%pythoncode %{
Geometries = Geometry()
%}

