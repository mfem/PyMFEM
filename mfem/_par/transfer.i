%module(package="mfem._par") transfer

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
%}

%init %{
import_array();
%}

%include "exception.i"
%import "operators.i"
%import "fespace.i"
%import "pfespace.i"
%include "../common/exception.i"

%import  "../common/pyoperator.hpp"
%include "fem/transfer.hpp"
