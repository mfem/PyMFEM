//
// Copyright (c) 2020-2026, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._par") ode
%{
#include  "mfem.hpp"
#include "linalg/ode.hpp"
#include "../common/pyoperator.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
%}

%init %{
import_array1(-1);
%}
%include "exception.i"
%import "vector.i"
%import "array.i"
%import "operators.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%typemap(in) double &t (double temp){
  temp = PyFloat_AsDouble($input);
  $1 = &temp;
 }
%typemap(in) double &dt (double dtemp){
  dtemp = PyFloat_AsDouble($input);
  $1 = &dtemp;
}
%typemap(argout) double &t {
  %append_output(PyFloat_FromDouble(*$1));
}
%typemap(argout) double &dt {
  %append_output(PyFloat_FromDouble(*$1));
 }

// Suppress the original wrappers using the protected alias.
// Then, inject PyMFEM version using publicly enum.
%ignore SupportsImplicitVariableType;
%ignore mfem::ODESolver::SetImplicitVariableType;

%rename(SupportsImplicitVariableType)
    mfem::ODESolver::PySupportsImplicitVariableType;
%rename(SetImplicitVariableType)
    mfem::ODESolver::PySetImplicitVariableType;

%include "linalg/ode.hpp"

%extend mfem::ODESolver {
    bool PySupportsImplicitVariableType(
        mfem::TimeDependentOperator::ImplicitVariableType var) const {
        return $self->SupportsImplicitVariableType(var);
    }

    void PySetImplicitVariableType(
        mfem::TimeDependentOperator::ImplicitVariableType var) {
        $self->SetImplicitVariableType(var);
    }
}


