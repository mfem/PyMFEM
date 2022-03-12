%module (package="mfem._par") auxiliary

%{
#include <iostream>  
#include "mfem.hpp"
#include "pyoperator.hpp"      
#include "numpy/arrayobject.h"    
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


