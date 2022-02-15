/*

   multigrid.i

*/
%module(package="mfem._par") multigrid
%feature("autodoc", "1");
%{
#include "linalg/operator.hpp"
#include "linalg/handle.hpp"
#include "fem/multigrid.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"    
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

LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Operator*>&,
			   mfem::Operator*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Solver*>&,
			   mfem::Solver*)
LIST_TO_MFEMOBJ_BOOLARRAY_IN(const mfem::Array<bool>& )
//%feature("director") mfem::PyGeometricMultigrid;

%pythonprepend mfem::PyGeometricMultigrid::AppendBilinearForm %{
   if not hasattr(self, "_forms"): self._forms = []
   self._forms.append(form)
   form.thisown = 0
%}
%pythonprepend mfem::PyGeometricMultigrid::AppendEssentialTDofs %{
   if not hasattr(self, "_esss"): self._esss = []
   self._esss.append(ess)	    
   ess.thisown = 0
%}
%feature("shadow") mfem::PyGeometricMultigrid::_pybfs %{
  @property						     
  def bfs(self):
     return self._forms
 %}       
%feature("shadow") mfem::PyGeometricMultigrid::_pyess %{       
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
  } /* end of namespace */ 
%}

  

