%module(package="mfem._par") vertex

%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}

%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%include "mesh/vertex.hpp"

