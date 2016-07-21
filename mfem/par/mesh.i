%module mesh
%{
#define MFEM_USE_MPI  
#include "mesh/mesh_headers.hpp"
#include "fem/fem.hpp"
#include "general/array.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
mfem::Mesh * MeshFromFile(const char *mesh_file, int generate_edges, int refine,
		      bool fix_orientation = true);
// void mfem:PrintToFile(const char *mesh_file,  const int precision) const;
#include "numpy/arrayobject.h"
#include "pycoefficient.hpp" 
%}
%init %{
import_array();
%}
%import "matrix.i"
%import "array.i"
%import "ncmesh.i"
%import "vector.i"
%import "element.i"
%import "vertex.i"
%import "mesh/mesquite.hpp"
%import "densemat.i"
%import "sparsemat.i"
%import "eltrans.i"
%import "intrules.i"
%feature("notabstract") VectorFunctionCoefficient;
%feature("notabstract") VectorConstantCoefficient;
%import "coefficient.i"
%import "fe.i"

//  conversion of Int (can handle numpy int)
%typemap(in) int {
  PyArray_PyIntAsInt($input);  
  $1 = PyInt_AsLong($input);
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) int {
  if (PyArray_PyIntAsInt($input)   != -1){
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%feature("shadow") mfem::Mesh::GetBdrElementVertices %{
def GetBdrElementVertices(self, i):
    from  .array import intArray
    ivert = intArray()
    _mesh.Mesh_GetBdrElementVertices(self, i, ivert)
    return ivert.ToList()
%}

%feature("shadow") mfem::Mesh::GetElementVertices %{
def GetElementVertices(self, i):
    from  .array import intArray
    ivert = intArray()
    _mesh.Mesh_GetElementVertices(self, i, ivert)
    return ivert.ToList()
%}

%feature("shadow") mfem::Mesh::GetElementEdges %{
def GetElementVEdges(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    _mesh.Mesh_GetElementEdges(self, i, ia, ib)
    return ia.ToList(), ib.ToList()      
%} 

%feature("shadow") mfem::Mesh::GetBdrElementEdges %{
def GetBdrElementEdges(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    _mesh.Mesh_GetBdrElementEdges(self, i, ia, ib)
    return ia.ToList(), ib.ToList()
%} 

%feature("shadow") mfem::Mesh::GetFaceEdges %{
def GetFaceEdges(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    _mesh.Mesh_GetFaceEdges(self, i, ia, ib)
    return ia.ToList(), ib.ToList()
%}

%feature("shadow") mfem::Mesh::GetEdgeVertices %{
def GetEdgeVertices(self, i):
    from  .array import intArray
    ia = intArray()
    _mesh.Mesh_GetEdgeVertices(self, i, ia)
    return ia.ToList()
%}

%feature("shadow") mfem::Mesh::GetFaceVertices %{
def GetFaceVertices(self, i):
    from  .array import intArray
    ia = intArray()
    _mesh.Mesh_GetFaceVertices(self, i, ia)
    return ia.ToList()
%}

%feature("shadow") mfem::Mesh::GetElementFaces %{
def GetElementFaces(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    _mesh.Mesh_GetElementFaces(self, i, ia, ib)
    return ia.ToList(), ib.ToList()
%}


%immutable attributes;
%immutable bdr_attributes;
%ignore MesquiteSmooth;
%include "mesh/mesh.hpp"
%mutable;
/*
%inline %{
#ifdef MFEM_USE_MPI
// auxiliary function for qsort
static int mfem_less(const void *x, const void *y)
{
   if (*(int*)x < *(int*)y)
   {
      return 1;
   }
   if (*(int*)x > *(int*)y)
   {
      return -1;
   }
   return 0;
}
// METIS 4 prototypes
typedef int idxtype;
extern "C" {  
   void METIS_PartGraphRecursive(int*, idxtype*, idxtype*, idxtype*, idxtype*,
                                 int*, int*, int*, int*, int*, idxtype*);
   void METIS_PartGraphKway(int*, idxtype*, idxtype*, idxtype*, idxtype*,
                            int*, int*, int*, int*, int*, idxtype*);
   void METIS_PartGraphVKway(int*, idxtype*, idxtype*, idxtype*, idxtype*,
                             int*, int*, int*, int*, int*, idxtype*);
}
#endif

 %}
*/
namespace mfem{
%extend Mesh{
   Mesh(const char *mesh_file, int generate_edges, int refine,
        bool fix_orientation = true){

        mfem::Mesh *mesh;
        std::ifstream imesh(mesh_file);
        if (!imesh)
        {
	  std::cerr << "\nCan not open mesh file: " << mesh_file << '\n' << std::endl;
   	  return NULL;
        }
	mesh = new mfem::Mesh(imesh, generate_edges, refine, fix_orientation);
	return mesh;
   }
   void PrintToFile(const char *mesh_file, const int precision) const
   {
	std::ofstream mesh_ofs(mesh_file);	
        mesh_ofs.precision(precision);
        self->Print(mesh_ofs);	
   }
   PyObject* GetVertexArray(int i) const
   {
     int L = self->Dimension();
     int n;
     const double *v = self->GetVertex(i);
     npy_intp dims[] = {L};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     double *x    = (double *)PyArray_DATA(array);
     for (n = 0; n < L; n++) {
        x[n] = v[n];
     }
     return array;
   }
   PyObject* GetBdrElementFace(int i) const
   {
     int a;
     int b;
     PyObject *o;
     
     if (i >= self->GetNBE()){
        return Py_BuildValue("");
     }
     self->GetBdrElementFace(i, &a, &b);
     o = Py_BuildValue("(ii)", a, b);
     return o;
   }
   PyObject* GetBdrAttributeArray() const
   {
     int i;
     npy_intp dims[] = {self->GetNBE()};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_INT);
     int *x    = (int *)PyArray_DATA(array);
     for (i = 0; i < self->GetNBE() ; i++){
       x[i] = self->GetBdrElement(i)->GetAttribute();
     }
     return array;
   }   
   PyObject* GetBdrArray(int idx) const
   {
     int i;
     int c = 0;
     for (i = 0; i < self->GetNBE() ; i++){
       if (self->GetBdrElement(i)->GetAttribute() == idx){c++;}
     }
     npy_intp dims[] = {c};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_INT);
     int *x    = (int *)PyArray_DATA(array);
     c = 0;     
     for (i = 0; i < self -> GetNBE() ; i++){
       if (self->GetBdrElement(i)->GetAttribute() == idx){
	 x[c] = (int)i;
         c++;
       }
     }
     return array;
   }   
  };   
}






