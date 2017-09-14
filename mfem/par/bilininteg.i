%module(directors="1")  bilininteg
%{
#include "fem/gridfunc.hpp"  
#include "fem/linearform.hpp"
#include "fem/bilininteg.hpp"
#include "pycoefficient.hpp"
#include "numpy/arrayobject.h"        
%}

%init %{
import_array();
%}

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
 //%template(IntegrationPointArray) mfem::Array<mfem::IntegrationPoint>;

%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
}

%include "../common/bilininteg_ext.i"
%include "fem/bilininteg.hpp"
