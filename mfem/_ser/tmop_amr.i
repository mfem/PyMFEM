%module(package="mfem._ser", directors="1") tmop_amr
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"  
%}

%include "exception.i"
%import "../common/exception_director.i"

%import tmop.i
%import nonlinearform.i

%include "fem/tmop_amr.hpp"
