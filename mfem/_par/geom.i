%module(package="mfem._par") geom
%{
#include "mfem.hpp"  
#include "fem/geom.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"  
%}

%init %{
import_array();
%}

%include "exception.i"
%import "intrules.i"
%import "densemat.i"

%immutable RefPts;
%immutable GlobGeometryRefiner;
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
  
   %template(GeometryTypeArray) Array<Geometry::Type>;
}

