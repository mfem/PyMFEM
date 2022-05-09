%module(package="mfem._par", directors="1")  bilininteg
%{
#include "mfem.hpp"
#include "pyoperator.hpp"    
#include "../common/pycoefficient.hpp"
#include "numpy/arrayobject.h"
  //using namespace mfem;
%}

%init %{
import_array();
%}

//%include "../common/cpointers.i"
//%import "cpointers.i"
%include "exception.i"

%import "globals.i"
%import "array.i"
%import "coefficient.i"
%import "matrix.i"
%import "vector.i"
%import "gridfunc.i"
%import "fespace.i"
%import "fe_coll.i"
%import "intrules.i"
%import "densemat.i"
%import "sparsemat.i"
%import "lininteg.i"
%import "eltrans.i"
%import "linearform.i"
%import "fe.i"
%import "nonlininteg.i"
%include "../common/exception_director.i"
 //%template(IntegrationPointArray) mfem::Array<mfem::IntegrationPoint>;

%feature("director") mfem::BilinearFormIntegrator;

%ignore  mfem::MassIntegrator::SetupPA;

%include "../common/bilininteg_ext.i"
%include "fem/bilininteg.hpp"
