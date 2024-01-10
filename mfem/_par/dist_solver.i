%module(package="mfem._par") dist_solver
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/common/dist_solver.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
%}

%init %{
import_array();
%}

%inline %{
#include "miniapps/common/dist_solver.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "pgridfunc.i"
%import "pmesh.i"
%import "solvers.i"

%include "miniapps/common/dist_solver.hpp"

