
%module mesh
%{
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
#include "iostream_typemap.hpp"         
#include "numpy/arrayobject.h"
#include "pycoefficient.hpp" 
%}

%init %{
import_array();
%}

%include "../common/cpointers.i"
%import "matrix.i"
%import "array.i"
%import "ncmesh.i"
%import "vector.i"
%import "gridfunc.i"
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
%import "ostream_typemap.i"
//

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
// this prevent automatic conversion from int to double so
// that it select collect overloaded method....
%typemap(typecheck,precedence=SWIG_TYPECHECK_DOUBLE) double {
  if (PyFloat_Check($input)){
    $1 = 1;
  } else {
    $1 = 0;
  }
}
// ignore these constructors, since in python element::type is given by 
// string (see extend section below).
// %ignore does not work well !?
//%ignore mfem::Mesh(int nx, int ny, int nz, mfem::Element::Type type,
//		   int generate_edges = 0, double sx = 1.0, double sy = 1.0,
//		   double sz = 1.0);
//%ignore mfem::Mesh(int nx, int ny, mfem::Element::Type type,
//                   int generate_edges = 0,
//		     double sx = 1.0, double sy = 1.0);
%typemap(typecheck) (int nx, int ny, int nz, mfem::Element::Type type) {
  $1 = 0; // ignore this pattern
}
%typemap(typecheck) (int nx, int ny, mfem::Element::Type type) {
  $1 = 0; // ignore this pattern
}

// to give vertex array as list
%typemap(in) (const double *){
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int l = PyList_Size($input);
  $1 = (double *) malloc((l)*sizeof(double));
  for (i = 0; i < l; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (PyInt_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
    } else if (PyFloat_Check(s)) {
        $1[i] = (double)PyFloat_AsDouble(s);
    } else {
        free($1);      
        PyErr_SetString(PyExc_ValueError, "List items must be integer/float");
        return NULL;
    }
  }
}
%typemap(typecheck) (const double *) {
   $1 = PyList_Check($input) ? 1 : 0;
}

// to give index array as list
%typemap(in) (const int *vi){
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int l = PyList_Size($input);
  $1 = (int *) malloc((l)*sizeof(int));
  for (i = 0; i < l; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyInt_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be integer");
        return NULL;
    }
    $1[i] = (int)PyInt_AsLong(s);
  }
}
%typemap(typecheck) (const int *vi) {
   $1 = PyList_Check($input) ? 1 : 0;
}

// SwapNodes


%typemap(in) mfem::GridFunction *&nodes (mfem::GridFunction *Pnodes){
int res2 = 0;
res2 = SWIG_ConvertPtr($input, (void **) &Pnodes, $descriptor(mfem::GridFunction *), 0);
if (!SWIG_IsOK(res2)){
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "Mesh_SwapNodes" "', argument " "2"" of type '" "*mfem::GridFunction""'");      
 }
 $1 = &Pnodes;
 }
 
%typemap(in) int &own_nodes_ (int own_nodes){
  own_nodes = (int)PyInt_AsLong($input);
  $1 = &own_nodes;
} 
%typemap(argout) (mfem::GridFunction *&nodes){
  Py_XDECREF($result);
  $result = PyList_New(0);
  %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*arg2), $descriptor(mfem::GridFunction *), 0 |  0 ));
 }
%typemap(argout) int &own_nodes_{
  %append_output(PyLong_FromLong((long)*$1));  
}


// default number is -1, which conflict with error code of PyArray_PyIntAsInt...
%typemap(typecheck) (int nonconforming = -1) {
   $1 = PyInt_Check($input) ? 1 : 0;
}


%feature("shadow") mfem::Mesh::GetBdrElementVertices %{
def GetBdrElementVertices(self, i):
    from  .array import intArray
    ivert = intArray()
    _mesh.Mesh_GetBdrElementVertices(self, i, ivert)
    return ivert.ToList()
%}

%feature("shadow") mfem::Mesh::GetBdrElementAdjacentElement %{
def GetBdrElementAdjacentElement(self, bdr_el):
    el = intp()
    info = intp()  
    _mesh.Mesh_GetBdrElementAdjacentElement(self, bdr_el, el, info)
    return el.value(), info.value()
%}

%feature("shadow") mfem::Mesh::GetElementVertices %{
def GetElementVertices(self, i):
    from  .array import intArray
    ivert = intArray()
    _mesh.Mesh_GetElementVertices(self, i, ivert)
    return ivert.ToList()
%}

%feature("shadow") mfem::Mesh::GetElementEdges %{
def GetElementEdges(self, i):
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

%feature("shadow") mfem::Mesh::GetBoundingBox %{
def GetBoundingBox(self, ref = 2):
    from  .vector import Vector
    min = Vector()
    max = Vector()      
    _mesh.Mesh_GetBoundingBox(self, min, max, ref)      
    return min.GetDataArray().copy(), max.GetDataArray().copy()
%}
%feature("shadow") mfem::Mesh::GetFaceElements %{
def GetFaceElements(self, Face):
    Elem1 = intp()
    Elem2 = intp()  
    val = _mesh.Mesh_GetFaceElements(self, Face, Elem1, Elem2)
    return Elem1.value(), Elem2.value()
%}

%immutable attributes;
%immutable bdr_attributes;
%ignore MesquiteSmooth;

%newobject mfem::Mesh::GetFaceToElementTable;
%newobject mfem::Mesh::GetVertexToElementTable;
%include "mesh/mesh.hpp"
%mutable;

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
   Mesh(int nx, int ny, int nz, const char *type, int generate_edges = 0,
        double sx = 1.0, double sy = 1.0, double sz = 1.0){
     mfem::Mesh *mesh;     
     if (std::strcmp(type, "POINT")) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::POINT,
			       generate_edges, sx, sy, sz);
     }
     else if (std::strcmp(type, "SEGMENT")) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::SEGMENT,
			       generate_edges, sx, sy, sz);
	 
     }
     else if (std::strcmp(type, "TRIANGLE")) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::TRIANGLE,
			       generate_edges, sx, sy, sz);
	 
     }
     else if (std::strcmp(type, "QUADRILATERAL")) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::QUADRILATERAL,
			       generate_edges, sx, sy, sz);
	 
     }	 
     else if (std::strcmp(type, "TETRAHEDRON")) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::TETRAHEDRON,
			       generate_edges, sx, sy, sz);
	 
     }	 
     else if (std::strcmp(type, "HEXAHEDRON")) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::HEXAHEDRON,
			       generate_edges, sx, sy, sz);
	 
     }	 
     else {
         return NULL;
     }
     return mesh;       
   }
   Mesh(int nx, int ny,  const char *type, int generate_edges = 0,
        double sx = 1.0, double sy = 1.0){
     mfem::Mesh *mesh;     
     if (std::strcmp(type, "POINT")) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::POINT,
			       generate_edges, sx, sy);
     }
     else if (std::strcmp(type, "SEGMENT")) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::SEGMENT,
			       generate_edges, sx, sy);
	 
     }
     else if (std::strcmp(type, "TRIANGLE")) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::TRIANGLE,
			       generate_edges, sx, sy);
	 
     }
     else if (std::strcmp(type, "QUADRILATERAL")) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::QUADRILATERAL,
			       generate_edges, sx, sy);
	 
     }	 
     else if (std::strcmp(type, "TETRAHEDRON")) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::TETRAHEDRON,
			       generate_edges, sx, sy);
	 
     }	 
     else if (std::strcmp(type, "HEXAHEDRON")) {
	 mesh = new mfem::Mesh(nx, ny,  mfem::Element::HEXAHEDRON,
			       generate_edges, sx, sy);
	 
     }	 
     else {
         return NULL;
     }
     return mesh;       
   }
   void PrintToFile(const char *mesh_file, const int precision) const
   {
	std::ofstream mesh_ofs(mesh_file);	
        mesh_ofs.precision(precision);
        self->Print(mesh_ofs);	
   }
   PyObject* GetAttributeArray() const
   {
     int i;
     npy_intp dims[] = {self->GetNE()};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_INT);
     int *x    = (int *)PyArray_DATA(array);
     for (i = 0; i < self->GetNE() ; i++){
       x[i] = (int)(self->GetElement(i)->GetAttribute());
     }
     return array;
   }   

   PyObject* GetVertexArray(int i) const
   {
     int L = self->SpaceDimension();     
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
       x[i] = (int)(self->GetBdrElement(i)->GetAttribute());
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
   PyObject* GetDomainArray(int idx) const
   {

     int i;
     int c = 0;     
     for (i = 0; i < self->GetNE() ; i++){
       if (self->GetElement(i)->GetAttribute() == idx){c++;}
     }
     npy_intp dims[] = {c};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_INT);
     int *x    = (int *)PyArray_DATA(array);
     c = 0;
     for (i = 0; i < self -> GetNE() ; i++){
       if (self->GetElement(i)->GetAttribute() == idx){
	 x[c] = (int)i;
         c++;
       }
     }
     return array;
   }
  };
}






