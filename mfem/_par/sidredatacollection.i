%module(package="mfem._par")  sidredatacollection
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
%import "datacollection.i"

%include "fem/sidredatacollection.hpp"
