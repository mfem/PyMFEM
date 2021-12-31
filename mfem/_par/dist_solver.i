%module(package="mfem._par") dist_solver
%{
#include "mfem.hpp"      
#include "miniapps/shifted/dist_solver.hpp"
#include "pyoperator.hpp"    
#include "../common/pycoefficient.hpp"  
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%inline %{
#include "miniapps/shifted/dist_solver.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "pgridfunc.i"
%import "pmesh.i"

%include "miniapps/shifted/dist_solver.hpp"

