/*

   multigrid.i

*/
%module(package="mfem._ser", director=1) multigrid
%feature("autodoc", "1");
%{
#include "linalg/operator.hpp"
#include "linalg/handle.hpp"
#include "fem/multigrid.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
%}
%init %{
import_array();
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "bilinearform.i"
%import "fespacehierarchy.i"
%import "../common/exception_director.i"
%import "../common/object_array_typemap.i"

ObjectArrayInput(mfem::Solver *);
ObjectArrayInput(mfem::Operator *);
BoolArrayInput(bool);

%feature("director") mfem::PyGeometricMultigrid;

%pythonprepend PyGeometricMultigrid::AppendBilinearForm %{
   if not hasattr(self, "_forms"): self._forms = []
   self._forms.append(form)
   form.thisown = 0
%}
%pythonprepend PyGeometricMultigrid::AppendEssentialTDofs %{
   if not hasattr(self, "_esss"): self._esss = []
   self._esss.append(ess)	    
   ess.thisown = 0
%}
%feature("shadow") PyGeometricMultigrid::_pybfs %{
  @property						     
  def bfs(self):
     return self._forms
 %}       
%feature("shadow") PyGeometricMultigrid::_pyess %{       
  @property						     
  def essentialTrueDofs(self):
     return self._esss
%}
%include "fem/multigrid.hpp"

%inline %{
  namespace mfem{
class PyGeometricMultigrid : public GeometricMultigrid
{
public:
 PyGeometricMultigrid(const FiniteElementSpaceHierarchy& fespaces_) 
   : GeometricMultigrid(fespaces_){}
  
  void AppendBilinearForm(BilinearForm *form){
    bfs.Append(form);
  }
  void AppendEssentialTDofs(Array<int> *ess){
      essentialTrueDofs.Append(ess);
  }
  void _pybfs(void){}
  void _pyess(void){}  
};
class CustomTransfer : public Operator
{
private:
  Operator *f_opr;
  SparseMatrix &M1;
  SparseMatrix &M2;
  
public:
 CustomTransfer(Operator *_f_opr, SparseMatrix &_M1, SparseMatrix &_M2)
   : Operator(_f_opr -> Height(), _f_opr -> Width()), f_opr(_f_opr), M1(_M1), M2(_M2){}
 
 virtual void Mult(const Vector &x, Vector &y) const {
   std::cout << "Forward \n";
   f_opr -> Mult(x, y);
 }
 virtual void MultTranspose(const Vector &x, Vector &y) const {
   std::cout << "Backward \n";
   Vector *B2 = new Vector(x.Size());
   Vector *B1 = new Vector(y.Size());
   
   M2.Mult(x, *B2);
   f_opr -> MultTranspose(*B2, *B1);

   GSSmoother prec(M1);
   PCG(M1, prec, *B1, y, -1, 2000, 1e-24, 0.0);

   delete B1;
   delete B2;
 }
};
  } /* end of namespace */
%}

  

