/*

   multigrid.i

*/
%module(package="mfem._par", director=1) multigrid
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
class PyGeometricMultigrid : public mfem::GeometricMultigrid
{
public:
 PyGeometricMultigrid(const mfem::FiniteElementSpaceHierarchy& fespaces_) 
   : mfem::GeometricMultigrid(fespaces_){}
  
  void AppendBilinearForm(mfem::BilinearForm *form){
    bfs.Append(form);
  }
  void AppendEssentialTDofs(mfem::Array<int> *ess){
      essentialTrueDofs.Append(ess);
  }
  void _pybfs(void){}
  void _pyess(void){}  
};
%}

  

