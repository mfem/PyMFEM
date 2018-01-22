%module(directors="1") operators
%{
#include "iostream_typemap.hpp"        
#include "linalg/operator.hpp"
#include "pyoperator.hpp"
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%import "vector.i"
%import "array.i"
%import "ostream_typemap.i"

%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
    //catch (...){
    //  SWIG_fail;
    //}
    //    catch (Swig::DirectorMethodException &e) { SWIG_fail; }
    //    catch (std::exception &e) { SWIG_fail; }    
}
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

%inline %{
void mfem::PyOperatorBase::Mult(const mfem::Vector &x, mfem::Vector &y) const
  {
    y = _EvalMult(x);
  }
void mfem::PyTimeDependentOperatorBase::Mult(const mfem::Vector &x, mfem::Vector &y) const
  {
    y = _EvalMult(x);
  }
%}

%feature("director") mfem::PyTimeDependentOperatorBase;
%feature("director") mfem::PyOperatorBase;
//%feature("noabstract") mfem::Operator;
//%feature("noabstract") mfem::TimeDependentOperator;
%feature("director") mfem::TimeDependentOperator;
%feature("director") mfem::Operator;

%include "linalg/operator.hpp"

%include "pyoperator.hpp"
%pythoncode %{
class PyOperator(PyOperatorBase):
   def __init__(self, *args):
       PyOperatorBase.__init__(self, *args)
   def _EvalMult(self, x):
       return self.EvalMult(x.GetDataArray())
   def EvalMult(self, x):
       raise NotImplementedError('you must specify this method')

class PyTimeDependentOperator(PyTimeDependentOperatorBase):
   def __init__(self, *args):  
       PyTimeDependentOperatorBase.__init__(self, *args)
   def _EvalMult(self, x):
       return self.EvalMult(x.GetDataArray())
   def EvalMult(self, x):
       raise NotImplementedError('you must specify this method')
			 
%}
  
