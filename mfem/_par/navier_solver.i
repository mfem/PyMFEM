%module(package="mfem._par") navier_solver
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/common/dist_solver.hpp"
// #include "miniapps/navier/navier_solver.hpp"
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
// #include "miniapps/navier/navier_solver.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "pgridfunc.i"
%import "pmesh.i"
%import "solvers.i"

%include "miniapps/common/dist_solver.hpp"
// %include "miniapps/navier/navier_solver.hpp"



// %module(package="mfem._par") navier_solver
// %{
// #include "mfem.hpp"
// #include "numpy/arrayobject.h"
// #include "miniapps/navier/navier_solver.hpp"
// #include "../common/pyoperator.hpp"
// #include "../common/pysolvers.hpp"
// #include "../common/pycoefficient.hpp"
// #include "../common/pyintrules.hpp"
// #include "../common/pybilininteg.hpp"
// #include "../common/pynonlininteg.hpp"
// %}

// %init %{
// import_array();
// %}

// %inline %{
// #include "miniapps/navier/navier_solver.cpp"
// %}


// %include "exception.i"
// %import "element.i"
// %import "../common/exception.i"

// %import "coefficient.i"
// %import "pgridfunc.i"
// %import "pmesh.i"
// %import "solvers.i"

// %include "miniapps/common/navier_solver.hpp"

