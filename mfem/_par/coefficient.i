/*
   coefficient.i
   SWIG interface file for coefficient.hpp

   Features
   1) function callback for VectorFunctionCoefficent


   Date: 2016. 2. 18
   Author: S. Shiraiwa (MIT)
 */
%module(package="mfem._par", directors="1")  coefficient
%feature("autodoc", "1");
/*%module  coefficient*/
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "mfem.hpp"  
#include "pyoperator.hpp"    
#include "../common/pycoefficient.hpp"
#include "numpy/arrayobject.h"
%}

%include "../common/mfem_config.i"
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%include "exception.i"

%import "globals.i"
%import "array.i"
%import "matrix.i"
%import "symmat.i"
%import "intrules.i"
%import "sparsemat.i"
%import "densemat.i"
%import "vector.i"
%import "eltrans.i"
%import "../common/exception_director.i"

%ignore Function;
//%ignore DeltaCoefficient;
%feature("notabstract") VectorFunctionCoefficient;
%feature("notabstract") VectorConstantCoefficient;
%feature("notabstract") VectorDeltaCoefficient;
%feature("notabstract") MatrixFunctionCoefficient;
%feature("notabstract") MatrixConstantCoefficient;
%feature("notabstract") mfem::CurlGridFunctionCoefficient;
%feature("notabstract") mfem::SymmetricMatrixConstantCoefficient;
%feature("notabstract") mfem::VectorQuadratureFunctionCoefficient;

namespace mfem {
%pythonprepend MatrixConstantCoefficient::MatrixConstantCoefficient(const DenseMatrix &m) %{
   try:
      import numpy as np
      value = np.array(m, copy=False, dtype=float)
      can_np_array = True
   except:
      can_np_array = False

   if can_np_array:
      v = mfem._par.vector.Vector(np.transpose(value).flatten())
      m = mfem._par.densemat.DenseMatrix(v.GetData(), value.shape[0], value.shape[1])       
      self._value = (v,m)
   else:
      pass 
%}
%pythonprepend VectorConstantCoefficient::VectorConstantCoefficient(const Vector &v) %{
   try:
      import numpy as np
      value = np.array(v, copy=False, dtype=float).flatten()
      can_np_array = True
   except:
      can_np_array = False

   if can_np_array:
      v = mfem._par.vector.Vector(value)
      self._value = v
   else:
      pass 
%}
%pythonprepend DeltaCoefficient::SetWeight %{
    w.thisown=0 
%}
%pythonprepend VectorArrayCoefficient::Set %{ 
    c.thisown=0 
%}
%pythonprepend MatrixArrayCoefficient::Set %{ 
    c.thisown=0 
%}
%pythonappend VectorRestrictedCoefficient::VectorRestrictedCoefficient %{
    self._ref_to_vc = vc
%}
%pythonappend RestrictedCoefficient::RestrictedCoefficient %{
    self._ref_to_c = c_
%}
%pythonappend MatrixRestrictedCoefficient::MatrixRestrictedCoefficient %{
    self._ref_to_mc = mc
%}
}

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

%include "../common/typemap_macros.i"
LIST_TO_MFEMOBJ_POINTERARRAY_IN(mfem::IntegrationRule const *irs[],  mfem::IntegrationRule *, 0)

%include "../../headers/coefficient.hpp"
%include "../common/numba_coefficient.i"

%feature("director") mfem::VectorPyCoefficientBase;
%feature("director") mfem::PyCoefficientBase;
%feature("director") mfem::MatrixPyCoefficientBase;

%inline %{
double fake_func(const mfem::Vector &x)
{
  return 0.0;
}
void fake_func_vec(const mfem::Vector &x, mfem::Vector &Ht)
{
     Ht(0) = 0.0;
     Ht(1) = 0.0;
     Ht(2) = 0.0;
}
void fake_func_mat(const mfem::Vector &x, mfem::DenseMatrix &Kt)
{
  Kt(0,0) = 1.0;
  Kt(1,0) = 0.0;
  Kt(2,0) = 0.0;
  Kt(0,1) = 0.0;
  Kt(1,1) = 1.0;
  Kt(2,1) = 0.0;
  Kt(0,2) = 0.0;
  Kt(1,2) = 0.0;
  Kt(2,2) = 1.0;
}

namespace mfem{ 
double PyCoefficientBase::Eval(ElementTransformation &T,
                               const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   if (isTimeDependent)
   {
     return _EvalPyT(transip, GetTime());
   }
   else
   {
     return _EvalPy(transip);
   }
}
void VectorPyCoefficientBase::Eval(Vector &V, ElementTransformation &T,
                                     const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   V.SetSize(vdim);
   if (isTimeDependent)
   {
      _EvalPyT(transip, GetTime(),  V);          
   }
   else
   {
      _EvalPy(transip, V);
   }
}

void VectorPyCoefficientBase::Eval(DenseMatrix &M, ElementTransformation &T,
                                  const IntegrationRule &ir)

{
   Vector Mi;
   M.SetSize(vdim, ir.GetNPoints());
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      M.GetColumnReference(i, Mi);
      const IntegrationPoint &ip = ir.IntPoint(i);
      T.SetIntPoint(&ip);
      Eval(Mi, T, ip);
   }
}

void MatrixPyCoefficientBase::Eval(DenseMatrix &K, ElementTransformation &T,
                                     const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);
   K.SetSize(height, width);   
   if (isTimeDependent)
   {
      _EvalPyT(transip, GetTime(),  K);          
   }
   else
   {
      _EvalPy(transip, K);
   }
}

}  /* end of name space*/
%}

%include "../common/pycoefficient.hpp"

%pythoncode %{
class PyCoefficient(PyCoefficientBase):
   def __init__(self):
       PyCoefficientBase.__init__(self, 0)
   def _EvalPy(self, x):
       return self.EvalValue(x.GetDataArray())
   def EvalValue(self, x):
       return 0.0
  
class PyCoefficientT(PyCoefficientBase):
   def __init__(self):
       PyCoefficientBase.__init__(self, 1)
   def _EvalPyT(self, x, t):
       return self.EvalValue(x.GetDataArray(), t)
   def EvalValue(self, x, t):
       return 0.0
	 
class VectorPyCoefficient(VectorPyCoefficientBase):
   def __init__(self, dim):
       self.vdim = dim
       VectorPyCoefficientBase.__init__(self, dim, 0)
   def _EvalPy(self, x, V):
       v = self.EvalValue(x.GetDataArray())
       V.Assign(v)	 

   def _EvalPyT(self, x, t, V):
       v = self.EvalValue(x.GetDataArray())
       V.Assign(v)	 	 

   def EvalValue(self, x):
       return [0,0,0]
  
class VectorPyCoefficientT(VectorPyCoefficientBase):
   def __init__(self, dim):
       self.vdim = dim  
       VectorPyCoefficientBase.__init__(self, dim, 1)
   def _EvalPy(self, x, V):
       v = self.EvalValue(x.GetDataArray(), 0)
       V.Assign(v)

   def _EvalPyT(self, x, t, V):
       v = self.EvalValue(x.GetDataArray(), t)
       V.Assign(v)	 	 	 

   def EvalValue(self, x, t):
       return [0.0,0.0,0.0]

class MatrixPyCoefficient(MatrixPyCoefficientBase):
   def __init__(self, dim):
       self.vdim = dim
       MatrixPyCoefficientBase.__init__(self, dim, 0)
   def _EvalPy(self, x, K):
       k = self.EvalValue(x.GetDataArray())
       K.Assign(k)	 	 	 	 	 	 

   def EvalValue(self, x):
       return np.array([[0,0,0], [0,0,0], [0,0,0]])
  
class MatrixPyCoefficientT(MatrixPyCoefficientBase):
   def __init__(self, dim):
       self.vdim = dim  
       MatrixPyCoefficientBase.__init__(self, dim, 1)
   def _EvalPyT(self, x, t, K):
       k = self.EvalValue(x.GetDataArray(), t)
       K.Assign(k)	 	 	 	 	 	 

   def EvalValue(self, x, t):
       return np.array([[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]])
	 
%}

