%module(package="mfem._par")  sidredatacollection
%{
#include "mfem.hpp"  
#include "pyoperator.hpp"        
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"  
%}

%init %{
import_array();
%}

%include "exception.i"
%import "datacollection.i"

%include "fem/sidredatacollection.hpp"
