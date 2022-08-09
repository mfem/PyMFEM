%module(package="mfem._ser", directors="1") tmop_tools
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}

%include "exception.i"
%import "../common/exception_director.i"

%include "../common/typemap_macros.i"

%import tmop.i
%import bilinearform.i

%include "fem/tmop_tools.hpp"
