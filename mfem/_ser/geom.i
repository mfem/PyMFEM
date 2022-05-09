%module(package="mfem._ser") geom
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"      
#include "../common/io_stream.hpp"      
%}

%init %{
import_array();
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

namespace mfem{
   %ignore Array<Geometry::Type>::Union;
   %ignore Array<Geometry::Type>::Find;
   %ignore Array<Geometry::Type>::FindSorted;
   %ignore Array<Geometry::Type>::Sort;
   %ignore Array<Geometry::Type>::DeleteFirst;
   %ignore Array<Geometry::Type>::Unique;
   %ignore Array<Geometry::Type>::PartialSum;
   %ignore Array<Geometry::Type>::Sum;
   %ignore Array<Geometry::Type>::IsSorted;
   %ignore Array<Geometry::Type>::Save;
   %ignore Array<Geometry::Type>::Max;
   %ignore Array<Geometry::Type>::Min;
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

