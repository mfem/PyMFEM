%module (package="mfem._par") auxiliary

%{
#include <iostream>  
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"
%include "exception.i"
%import "../common/exception.i"

%import "solvers.i"
%import "pfespace.i"

%include "linalg/auxiliary.hpp"


