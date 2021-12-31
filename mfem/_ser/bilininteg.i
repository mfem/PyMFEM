%module(package="mfem._ser", directors="1")  bilininteg
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

/*
%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
    //catch (...){
    //  SWIG_fail;
    //}
    //    catch (Swig::DirectorMethodException &e) { SWIG_fail; }
    //    catch (std::exception &e) { SWIG_fail; }    
}
*/
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

 //%feature("director:except") {
 //    if ($error != NULL) {
 //        throw Swig::DirectorMethodException();
 //    }
 //}
%feature("director") mfem::BilinearFormIntegrator;

%include "../common/bilininteg_ext.i"

%ignore  mfem::MassIntegrator::SetupPA;

%include "fem/bilininteg.hpp"
