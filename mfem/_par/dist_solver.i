%module(package="mfem._par") dist_solver
%{
#include "config/config.hpp"      
#include "fem/fem.hpp"
#include "fem/coefficient.hpp"  
#include "miniapps/shifted/dist_solver.hpp"
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

