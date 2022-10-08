%module(package="mfem._ser", directors="0")  mesh
  
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include <vector>
#include "mfem.hpp"  
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/io_stream.hpp"
%}

%begin %{
#define PY_SSIZE_T_CLEAN
%}
%init %{
import_array();
%}

%include "exception.i"

%include "std_string.i"

%import "matrix.i"
%import "mem_manager.i"
%import "array.i"
%import "sort_pairs.i"
%import "ncmesh.i"
%import "vector.i"
%import "gridfunc.i"
%import "element.i"
%import "vertex.i"
%import "vtk.i"
%import "mesh/mesquite.hpp"
%import "densemat.i"
%import "sparsemat.i"
%import "eltrans.i"
%import "intrules.i"

%import "std_vectors.i"

%feature("notabstract") VectorFunctionCoefficient;
%feature("notabstract") VectorConstantCoefficient;
%import "coefficient.i"
%import "fe.i"
%import "../common/numpy_int_typemap.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

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

%import "../common/const_doubleptr_typemap.i"
CONST_DOUBLEPTR_IN(const double *)

%import "../common/const_intptr_typemap.i"
CONST_INTPTR_IN(const int *vi)


// to give index array as list/Array<int>/const int*

// SwapNodes (
//   it return new *GridFunction and own_nodes, also if nodes is NULL
//   it return None
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
  if (*arg$argnum){
     // return None if Nodes is NULL
     %append_output(Py_None);
  } else {
     %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(*arg$argnum), $descriptor(mfem::GridFunction *), 0 |  0 ));
  }
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
    $action(self, i, ivert)
    return ivert.ToList()
%}

%feature("shadow") mfem::Mesh::GetBdrElementAdjacentElement %{
def GetBdrElementAdjacentElement(self, bdr_el):
    from mfem.ser import intp
    el = intp()
    info = intp()  
    $action(self, bdr_el, el, info)
    return el.value(), info.value()
%}

%feature("shadow") mfem::Mesh::GetElementVertices %{
def GetElementVertices(self, i):
    from  .array import intArray
    ivert = intArray()
    $action(self, i, ivert)
    return ivert.ToList()
%}

%feature("shadow") mfem::Mesh::GetElementEdges %{
def GetElementEdges(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    $action(self, i, ia, ib)
    return ia.ToList(), ib.ToList()      
%} 

%feature("shadow") mfem::Mesh::GetBdrElementEdges %{
def GetBdrElementEdges(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    $action(self, i, ia, ib)
    return ia.ToList(), ib.ToList()
%} 

%feature("shadow") mfem::Mesh::GetFaceEdges %{
def GetFaceEdges(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    $action(self, i, ia, ib)
    return ia.ToList(), ib.ToList()
%}

%feature("shadow") mfem::Mesh::GetEdgeVertices %{
def GetEdgeVertices(self, i):
    from  .array import intArray
    ia = intArray()
    $action(self, i, ia)
    return ia.ToList()
%}

%feature("shadow") mfem::Mesh::GetFaceVertices %{
def GetFaceVertices(self, i):
    from  .array import intArray
    ia = intArray()
    $action(self, i, ia)
    return ia.ToList()
%}

%feature("shadow") mfem::Mesh::GetElementFaces %{
def GetElementFaces(self, i):
    from  .array import intArray
    ia = intArray()
    ib = intArray()      
    $action(self, i, ia, ib)
    return ia.ToList(), ib.ToList()
%}

%feature("shadow") mfem::Mesh::GetBoundingBox %{
def GetBoundingBox(self, ref = 2):
    from  .vector import Vector
    min = Vector()
    max = Vector()      
    $action(self, min, max, ref)      
    return min.GetDataArray().copy(), max.GetDataArray().copy()
%}

%feature("shadow") mfem::Mesh::GetFaceElements %{
def GetFaceElements(self, Face):
    from mfem.ser import intp
    Elem1 = intp()
    Elem2 = intp()  
    val = $action(self, Face, Elem1, Elem2)
    return Elem1.value(), Elem2.value()
%}
%feature("shadow") mfem::Mesh::GetElementTransformation %{
def GetElementTransformation(self, i):
    from mfem.ser import IsoparametricTransformation
    Tr = IsoparametricTransformation()
    $action(self, i, Tr)
    return Tr
%}
%feature("shadow") mfem::Mesh::GetBdrElementTransformation %{
def GetBdrElementTransformation(self, i):
    from mfem.ser import IsoparametricTransformation
    Tr = IsoparametricTransformation()
    $action(self, i, Tr)
    return Tr
%}
%feature("shadow") mfem::Mesh::GetFaceTransformation %{
def GetFaceTransformation(self, i):
    from mfem.ser import IsoparametricTransformation
    Tr = IsoparametricTransformation()
    $action(self, i, Tr)
    return Tr
%}
%feature("shadow") mfem::Mesh::GetEdgeTransformation %{
def GetEdgeTransformation(self, i):
    from mfem.ser import IsoparametricTransformation
    Tr = IsoparametricTransformation()
    $action(self, i, Tr)
    return Tr
%}
%feature("shadow") mfem::Mesh::GetFaceInfos %{
def GetFaceInfos(self, i):
    from mfem.ser import intp
    Elem1 = intp()
    Elem2 = intp()  
  
    $action(self, i, Elem1, Elem2)
    return Elem1.value(), Elem2.value()
%}
%feature("shadow") mfem::Mesh::FindPoints %{
def FindPoints(self, pp, warn=True, inv_trans=None):          
    r"""count, element_id, integration_points = FindPoints(points, warn=True, inv_trans=None)"""
    import numpy as np
    import mfem.ser as mfem      

    pp = np.array(pp, copy=False, dtype=float).transpose()      
    M = mfem.DenseMatrix(pp.shape[0], pp.shape[1])
    M.Assign(pp)
    elem_ids = mfem.intArray()
    int_points = mfem.IntegrationPointArray()
    count = $action(self, M, elem_ids, int_points, warn, inv_trans)
    elem_ids = elem_ids.ToList()
    return count, elem_ids, int_points
%}
%feature("shadow") mfem::Mesh::CartesianPartitioning %{
def CartesianPartitioning(self, nxyz, return_list=False):
    import mfem.ser as mfem
    import warnings      
    try:
        nxyz = list(nxyz)
        d = mfem.intArray(nxyz)
        dd = d.GetData()
    except BaseException:
        dd = nxyz
        warnings.warn("CartesianPartitioning argument should be iterable",
		      DeprecationWarning,)
    r = $action(self, dd)

    if not return_list:
        return r
    else:	 
        result = mfem.intArray()
        result.MakeRef(r, self.GetNE())
        result.MakeDataOwner()
        return result.ToList()
%}

%immutable attributes;
%immutable bdr_attributes;
%ignore MesquiteSmooth;

%newobject mfem::Mesh::GetFaceToElementTable;
%newobject mfem::Mesh::GetVertexToElementTable;

%include "../common/exception.i"
%include "mesh/mesh.hpp"

%mutable;

namespace mfem{
%extend Mesh{
    //Mesh(const char *mesh_file, int generate_edges, int refine,
    //    bool fix_orientation = true){
    //
    //    mfem::Mesh *mesh;
    //    std::ifstream imesh(mesh_file);
    //    if (!imesh)
    //    {
    //	  std::cerr << "\nCan not open mesh file: " << mesh_file << '\n' << std::endl;
    //  return NULL;
    //    }
    //	mesh = new mfem::Mesh(imesh, generate_edges, refine, fix_orientation);
    //	return mesh;
    //   }
   Mesh(int nx, int ny, int nz, const char *type, bool generate_edges = 0,
        double sx = 1.0, double sy = 1.0, double sz = 1.0,
	bool sfc_ordering = true){
     mfem::Mesh *mesh;     

     if (std::strcmp(type, "TETRAHEDRON") == 0) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::TETRAHEDRON,
			       generate_edges, sx, sy, sz);
	 
     }	 
     else if (std::strcmp(type, "HEXAHEDRON") == 0) {
	 mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::HEXAHEDRON,
			       generate_edges, sx, sy, sz);
	 
     }	 
     else {
         return NULL;
     }
     return mesh;       
   }
   Mesh(int nx, int ny,  const char *type, bool generate_edges = 0,
        double sx = 1.0, double sy = 1.0, bool sfc_ordering = true){
     mfem::Mesh *mesh;
     if (std::strcmp(type, "TRIANGLE") == 0) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::TRIANGLE,
			       generate_edges, sx, sy);
	 
     }
     else if (std::strcmp(type, "QUADRILATERAL") == 0) {
	 mesh = new mfem::Mesh(nx, ny, mfem::Element::QUADRILATERAL,
			       generate_edges, sx, sy);
	 
     }	 
     else {
         return NULL;
     }
     return mesh;       
   }

   static PyObject * MakeMerged(PyObject *tuple_or_list){
     PyObject *input = tuple_or_list;
     bool is_tuple = false;

     if (PyList_Check(input)) {
       is_tuple = false;
     }
     else if (PyTuple_Check(input)) {
       is_tuple = true;
     }
     else {
       PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
       return NULL;
     }

     int size = (is_tuple) ? PyTuple_Size(input) : PyList_Size(input);

     mfem::Mesh *mesh_array[size];
     mfem::Mesh *tmp_ptr;

     swig_type_info *ty = $descriptor(const mfem::Mesh *);

     for (int i = 0; i < size; i++) {
       PyObject *s = is_tuple ? PyTuple_GetItem(input, i) : PyList_GetItem(input,i);
       if (!SWIG_IsOK(SWIG_ConvertPtr(s, (void **)&tmp_ptr, ty, 0|0))) {
	 PyErr_SetString(PyExc_ValueError, "List items must be mfem::Mesh *");
         return NULL;
       } else {
	 mesh_array[i] = tmp_ptr;
       }
     }
     mfem::Mesh *mesh = new mfem::Mesh(mesh_array, size);
     int own = 1;
     return SWIG_NewPointerObj(SWIG_as_voidptr(mesh), ty, own);
   }

   void PrintToFile(const char *mesh_file, const int precision) const
   {
        std::cerr << "\nWarning Deprecated : Use Print(filename) insteead of SaveToFile \n";     
	std::ofstream mesh_ofs(mesh_file);	
        mesh_ofs.precision(precision);
        self->Print(mesh_ofs);	
   }

   PyObject* WriteToStream(PyObject* StringIO) const
   {
      PyObject* module = PyImport_ImportModule("io");
      if (!module){
   	 PyErr_SetString(PyExc_RuntimeError, "Can not load io module");
         return (PyObject *) NULL;
      }      
      PyObject* cls = PyObject_GetAttrString(module, "StringIO");
      if (!cls){
   	 PyErr_SetString(PyExc_RuntimeError, "Can not load StringIO");
         return (PyObject *) NULL;
      }      
      int check = PyObject_IsInstance(StringIO, cls);
      Py_DECREF(module);
      if (! check){
 	 PyErr_SetString(PyExc_TypeError, "First argument must be IOString");
         return (PyObject *) NULL;
      }
      std::ostringstream stream;
      self->Print(stream);      
      std::string str =  stream.str();
      const char* s = str.c_str();
      const int n = str.length();
      PyObject *ret = PyObject_CallMethod(StringIO, "write", "s#", s, static_cast<Py_ssize_t>(n));
      if (PyErr_Occurred()) {
         PyErr_SetString(PyExc_RuntimeError, "Error occured when writing IOString");
         return (PyObject *) NULL;
      }
      return ret;
   }
   
   PyObject* GetAttributeArray() const
   {
     int i;
     npy_intp dims[] = {self->GetNE()};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_INT);
     int *x    = (int *)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
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
     return  PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (double *)v);
   }

   PyObject* GetVertexArray() const
   {
     int L = self->SpaceDimension();
     int NV = self->GetNV();          
     const double *data;
     npy_intp dims[] = {L};

     PyObject *list = PyTuple_New(NV);
     for (int i = 0; i < NV; i++) {
          data = self->GetVertex(i);
          PyObject *array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (double *)data);
          PyTuple_SetItem(list, i, array);
     }
     return list;
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
     int *x    = (int *)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
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
     int *x    = (int *)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
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
     int *x    = (int *)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
     c = 0;
     for (i = 0; i < self -> GetNE() ; i++){
       if (self->GetElement(i)->GetAttribute() == idx){
	 x[c] = (int)i;
         c++;
       }
     }
     return array;
   }
   
   PyObject* GetElementCenterArray(int idx)
   {
     int i;
     mfem::Vector v;

     self->GetElementCenter(idx, v);

     npy_intp dims[] = {v.Size()};
     PyObject *array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
     double *x    = (double *)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
     for (i = 0; i < v.Size() ; i++){
	 x[i] = v[i];
     }
     return array;
   }

  double GetScaledJacobian(int i, int sd=2)
  {
    // compute scaled Jacobian
    //   i : element index
    //   sd: subdivision 
    //   https://github.com/mfem/mfem/pull/1835/files
    // 
    double attr = mfem::infinity();
    mfem::DenseMatrix J(self->Dimension());
    
    mfem::Geometry::Type geom = self->GetElementBaseGeometry(i);
    mfem::ElementTransformation *T = self->GetElementTransformation(i);

    mfem::RefinedGeometry *RefG = mfem::GlobGeometryRefiner.Refine(geom, sd, 1);
    mfem::IntegrationRule &ir = RefG->RefPts;

    // For each element, find the minimal scaled Jacobian in a
    // lattice of points with the given subdivision factor.

    for (int j = 0; j < ir.GetNPoints(); j++)
      {
	T->SetIntPoint(&ir.IntPoint(j));
	mfem::Geometries.JacToPerfJac(geom, T->Jacobian(), J);

        // Jacobian determinant
        double sJ = J.Det();

        for (int k = 0; k < J.Width(); k++)
	  {
	    mfem::Vector col;
	    J.GetColumnReference(k, col);
            // Scale by column norms
	    sJ /= col.Norml2();
	  }

	attr = fmin(sJ, attr);
      }
    return attr;
  }

  PyObject *IsElementOnPlaneArray(PyObject *aa, PyObject *bb, PyObject *cc, PyObject *dd){  
    /*
    return Boolean numpy array to indicate which element is on the plane
    defined by ax + by + cz + d = 0
    */
    double a = PyFloat_AsDouble(aa);    
    double b = PyFloat_AsDouble(bb);
    double c = PyFloat_AsDouble(cc);
    double d = PyFloat_AsDouble(dd);    
    /*
    return Boolean numpy array to indicate which element is on the plane
    defined by ax + by + cz + d = 0
    */
    mfem::Array<int> inodes;
    int nele = self -> GetNE();
    double *ptx;
    npy_intp dims[] = {nele};
    PyObject *array = PyArray_SimpleNew(1, dims, NPY_BOOL);
    if (self -> SpaceDimension() != 3  || self -> Dimension() != 3){
 	 PyErr_SetString(PyExc_TypeError, "dim and sdim must be 3");
         return (PyObject *) NULL;
    }
    bool *x    = (bool *)PyArray_DATA(reinterpret_cast<PyArrayObject *>(array));
    bool check, check2;
    
    for (int k = 0; k < self -> GetNE(); k++){
      self-> GetElementVertices(k, inodes);
      ptx = self -> GetVertex(inodes[0]);    
      check = (ptx[0]*a + ptx[1]*b + ptx[2]*c + d >=0);
      x[k] = false;
      for (int j=0; j < inodes.Size(); j++){
	ptx = self -> GetVertex(inodes[j]);
        check2 = (ptx[0]*a + ptx[1]*b + ptx[2]*c + d >=0);
        if (check != check2){
	  x[k] = true;
	  break;
	}
      }
    }
    return array;
  }
  };  // end of extend
}     // end of namespace

/*
virtual void PrintXG(std::ostream &out = mfem::out) const;
virtual void Print(std::ostream &out = mfem::out) const { Printer(out); }
void PrintVTK(std::ostream &out);
virtual void PrintInfo(std::ostream &out = mfem::out)
*/

OSTREAM_ADD_DEFAULT_FILE(Mesh, PrintInfo)
OSTREAM_ADD_DEFAULT_FILE(Mesh, Print)
OSTREAM_ADD_DEFAULT_FILE(Mesh, PrintXG)
OSTREAM_ADD_DEFAULT_FILE(Mesh, PrintVTK)



