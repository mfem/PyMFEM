%module(package="mfem._par") sets
%{
#include <mpi.h>
#include "general/sets.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%include "exception.i"
%include "../common/exception.i"

%import "array.i"
%import "table.i"
%include "general/sets.hpp"
