%module(package="mfem._ser") datacollection
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}
%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%import "globals.i"
%import "mesh.i"
%import "gridfunc.i"

%include "fem/datacollection.hpp"
