%module(directors="1")  lininteg
%{
#include "fem/lininteg.hpp"  
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "pycoefficient.hpp"
%}

 //%include "fem/coefficient.hpp"
%import "fe.i"
%import "vector.i"
%import "eltrans.i"
%import "intrules.i"
%import "coefficient.i"


%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
    //catch (...){
    //  SWIG_fail;
    //}
    //    catch (Swig::DirectorMethodException &e) { SWIG_fail; }
    //    catch (std::exception &e) { SWIG_fail; }    
}

namespace mfem { 
%pythonappend VectorFEBoundaryTangentLFIntegrator::VectorFEBoundaryTangentLFIntegrator %{
    self._coeff = QG
%}
%pythonappend VectorFEDomainLFIntegrator::VectorFEDomainLFIntegrator %{
   self._coeff = F
%}
}

%include "fem/lininteg.hpp"


