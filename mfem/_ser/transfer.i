%module(package="mfem._ser") transfer

%feature("autodoc", "1");

%{
#include  "mfem.hpp"
#include "mesh/vertex.hpp"
#include "numpy/arrayobject.h"
#include "fem/transfer.hpp"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp" 
%}

%init %{
import_array();
%}

%include "exception.i"
%import "operators.i"
%import "fespace.i"
%include "../common/exception.i"

%import  "pyoperator.hpp"
%include "fem/transfer.hpp"

