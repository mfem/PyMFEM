%module(package="mfem._ser") sort_pairs
%{
#include  "mfem.hpp"
#include "general/globals.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}
%include "exception.i"
%include "../common/typemap_macros.i"
%include "../common/exception.i"

%include "general/sort_pairs.hpp"

namespace mfem{
%template(intintintTriple) Triple<int, int, int>;
 }
