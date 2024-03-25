%module(package="mfem._par") navier_solver
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "miniapps/navier/navier_solver.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%init %{
import_array();
%}

%inline %{
#include "miniapps/navier/navier_solver.cpp"
%}


%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "coefficient.i"
%import "pgridfunc.i"
%import "pmesh.i"
%import "solvers.i"
%import "operators.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

OSTREAM_TYPEMAP(std::ostream&)


%typemap(in) double &time (double temp){
  temp = PyFloat_AsDouble($input);
  $1 = &temp;
 }

%typemap(argout) double &time {
  %append_output(PyFloat_FromDouble(*$1));
}

%ignore *::AddVelDirichletBC(VecFuncT *f, Array &);
%ignore *::AddPresDirichletBC(ScalarFuncT *f, Array &);
%ignore *::AddAccelTerm(VecFuncT *f, Array &);

%include "miniapps/navier/navier_solver.hpp"

