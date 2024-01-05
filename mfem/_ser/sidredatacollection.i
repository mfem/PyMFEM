%module(package="mfem._ser")  sidredatacollection
%{
#include "mfem.hpp"  
#include "../common/pyoperator.hpp"        
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"  
%}

%init %{
import_array();
%}

%include "exception.i"
%import "datacollection.i"

%include "fem/sidredatacollection.hpp"
