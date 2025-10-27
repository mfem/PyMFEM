//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
/*

   multigrid.i

*/
%module(package="mfem._ser") multigrid
%feature("autodoc", "1");
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
%}
%init %{
import_array1(-1);
%}

%include "exception.i"
%import "array.i"
%import "vector.i"
%import "bilinearform.i"
%import "fespacehierarchy.i"
%import "../common/exception_director.i"
//%import "../common/object_array_typemap.i"
%include "../common/typemap_macros.i"

LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Operator*>&,
			   mfem::Operator*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Solver*>&,
			   mfem::Solver*)
LIST_TO_MFEMOBJ_BOOLARRAY_IN(const mfem::Array<bool>& )

%pythonprepend mfem::PyGeometricMultigrid::AppendBilinearForm %{
  if not hasattr(self, "forms"):
      self.forms=[]
  self.forms.append(form)
  form.thisown=0
%}

%feature("shadow") mfem::PyGeometricMultigrid::_pybfs %{
  @property
  def bfs(self):
     return self._GetBilinearFormArray()
%}

%feature("shadow") mfem::PyGeometricMultigrid::_pyess %{
  @property
  def essentialTrueDofs(self):
     return self._GetEssentialTrueDofs()
%}

%include "fem/multigrid.hpp"

%inline %{
  namespace mfem{
class PyGeometricMultigrid : public GeometricMultigrid
{
public:
 PyGeometricMultigrid(const FiniteElementSpaceHierarchy& fespaces_)
   : GeometricMultigrid(fespaces_){}
 PyGeometricMultigrid(const FiniteElementSpaceHierarchy& fespaces_,
                      const Array<int> &ess_bdr_)
   : GeometricMultigrid(fespaces_, ess_bdr_){}

  void AppendBilinearForm(BilinearForm *form){
    bfs.Append(form);
  }

  Array<BilinearForm*> *_GetBilinearFormArray(){
    return &bfs;
  }

  Array<Array<int> *> *_GetEssentialTrueDofs(){
    return &essentialTrueDofs;
  }

  void _pybfs(void){}
  void _pyess(void){}
};
  } /* end of namespace */
%}
