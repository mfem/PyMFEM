%module fespace
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>

#include "iostream_typemap.hpp"            
  
#include "numpy/arrayobject.h"
#include "fem/fem.hpp"
#include "fem/fe_coll.hpp"
#include "fem/fespace.hpp"
#include "fem/eltrans.hpp"
#include "fem/coefficient.hpp"
#include "fem/intrules.hpp"  
#include "pyoperator.hpp"       
%}

%init %{
import_array();
%}

%import "array.i"
%import "vector.i"
%import "coefficient.i"
%import "matrix.i"
%import "mesh.i"
%import "intrules.i"
%import "fe.i"
%import "fe_coll.i"
%import "densemat.i"
%import "sparsemat.i"
%import "eltrans.i"
%import "lininteg.i"
%import "ostream_typemap.i"

// default number is -1, which conflict with error code of PyArray_PyIntAsInt...
%typemap(typecheck) (int ndofs = -1) {
   $1 = PyInt_Check($input) ? 1 : 0;
}
%typemap(typecheck) (int component = -1) {
   $1 = PyInt_Check($input) ? 1 : 0;
}


//VDoF accesser
%feature("shadow") mfem::FiniteElementSpace::GetBdrElementVDofs %{
def GetBdrElementVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetBdrElementVDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementVDofs %{
def GetElementVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetElementVDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetFaceVDofs %{
def GetFaceVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetFaceVDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeVDofs %{
def GetEdgeVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetEdgeVDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetVertexVDofs %{
def GetVertexVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetVertexVDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementInteriorVDofs %{
def GetElementInteriorVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetElementInteriorVDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeInteriorVDofs %{
def GetEdgeInteriorVDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetEdgeInteriorVDofs(self, i, vdofs)
    return vdofs.ToList()
%}

//DoF accesser
%feature("shadow") mfem::FiniteElementSpace::GetBdrElementDofs %{
def GetBdrElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetBdrElementDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementDofs %{
def GetElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetElementDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetFaceDofs %{
def GetFaceDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetFaceDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeDofs %{
def GetEdgeDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetEdgeDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetVertexDofs %{
def GetVertexDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetVertexDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetElementInteriorDofs %{
def GetElementInteriorDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetElementInteriorDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::FiniteElementSpace::GetEdgeInteriorDofs %{
def GetEdgeInteriorDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _fespace.FiniteElementSpace_GetEdgeInteriorDofs(self, i, vdofs)
    return vdofs.ToList()
%}

%pythonappend mfem::FiniteElementSpace::FiniteElementSpace%{

    '''
        FiniteElementSpace(Mesh *m, const FiniteElementCollection *f,
                      int vdim = 1, int ordering = Ordering::byNODES);
        keep reference to mesh and fec to prevent it from
        freeed from pytho garbage collector
    '''
    if hasattr(self, 'this'):
        self.mesh = args[0]
        self.fec = args[1]
      
%}  
%include "fem/fespace.hpp"

