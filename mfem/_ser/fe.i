%module(package="mfem._ser") fe

%{
#include <iostream>  
#include "mfem.hpp"
#include "pyoperator.hpp"      
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%immutable;
%ignore poly1d;
%mutable;

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "geom.i"
%import "intrules.i"
%import "densemat.i"
%import "sparsemat.i"
%import "fe_base.i"
%import "fe_fixed_order.i"
%import "fe_h1.i"
%import "fe_nd.i"
%import "fe_rt.i"
%import "fe_l2.i"
%import "fe_nurbs.i"
%import "fe_pos.i"
%import "fe_ser.i"
%import "../common/exception.i"

%ignore mfem::DofToQuad::FE;
%include "fem/fe.hpp"

