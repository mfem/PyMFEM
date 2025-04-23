%module(package="mfem._par") fe_base
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
%import "intrules.i"
%import "geom.i"
%import "doftrans.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

//forward declearation
%inline %{
namespace mfem{
  class FiniteElement;
}
%}

// this is to avoild extern breaks template
namespace mfem{
  %ignore poly1d;
}

%include "fem/fe/fe_base.hpp"
