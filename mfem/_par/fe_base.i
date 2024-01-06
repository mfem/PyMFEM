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

 //%ignore FE;
namespace mfem{
  class FiniteElement;
}

%include "fem/fe/fe_base.hpp"
