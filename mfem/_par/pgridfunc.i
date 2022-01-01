%module(package="mfem._par") pgridfunc

%feature("autodoc", "1");

%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>  
#include <mpi.h>
#include "mfem.hpp"
#include "pyoperator.hpp"  
#include "../common/io_stream.hpp"      
#include "../common/pycoefficient.hpp"  
#include "numpy/arrayobject.h"
%}

%include "../common/mfem_config.i"

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%init %{
import_array();
%}

%include "exception.i"
 //%include "../common/cpointers.i"
 //%import "cpointers.i"
%import "pfespace.i"
%import "gridfunc.i"
%import "hypre.i"
%import "pmesh.i"
%import "linearform.i"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%pointer_class(int, intp);

/*
%typemap(in) const mfem::IntegrationRule *irs[]{
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (mfem::IntegrationRule **) malloc((size)*sizeof(mfem::IntegrationRule *));
    for (i = 0; i < size; i++) {
       PyObject *o = PyList_GetItem($input,i);
       void *temp;       
       if (SWIG_ConvertPtr(o, &temp,
	   $descriptor(mfem::IntegrationRule *),SWIG_POINTER_EXCEPTION) == -1){
           return NULL;
       }
       $1[i] = reinterpret_cast<mfem::IntegrationRule *>(temp);       
     }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap(typecheck) const mfem::IntegrationRule *irs[]{
   $1 = PyList_Check($input) ? 1 : 0;
}
*/
//%rename(Assign) mfem::ParGridFunction::operator=;

%newobject mfem::ParGridFunction::ParallelAssemble;
%newobject mfem::ParGridFunction::ParallelAverage;
%newobject mfem::ParGridFunction::ParallelProject;
%newobject mfem::ParGridFunction::GetTrueDofs;


%pythonprepend mfem::ParGridFunction::ParGridFunction %{
from mfem._par.pmesh import ParMesh
from mfem._par.pfespace import ParFiniteElementSpace
from mfem._par.gridfunc import GridFunction
if (len(args) == 2 and isinstance(args[1], str) and
     isinstance(args[0], ParMesh)):
    g0 = GridFunction(args[0], args[1])
    fec = g0.OwnFEC()
    fes = g0.FESpace()
    pfes = ParFiniteElementSpace(args[0], fec, fes.GetVDim(),
                                      fes.GetOrdering())
    x = ParGridFunction(pfes, g0)
    x.thisown = 0
    pfes.thisown = 0
    g0.thisown = 0
    self.this = x.this
    return 
%}
   
%include "../common/typemap_macros.i"
LIST_TO_MFEMOBJ_POINTERARRAY_IN(mfem::IntegrationRule const *irs[],  mfem::IntegrationRule *, 0)

%include "fem/pgridfunc.hpp"

namespace mfem{
%extend ParGridFunction{
ParGridFunction(mfem::ParFiniteElementSpace *fes, const mfem::Vector &v, int offset){
   mfem::ParGridFunction *gf;   
   gf = new mfem::ParGridFunction(fes, v.GetData() + offset);
   return gf;
 }
  void Assign(const mfem::HypreParVector &v) {
    (* self) = v;
  }
  void Assign(const mfem::ParGridFunction &v) {
    (* self) = v;
  }
  void Assign(const double v) {
    (* self) = v;
  }
  void Assign(const mfem::Vector &v) {
    (* self) = v;
  }
  //void Assign(const mfem::GridFunction &v) {
  //  (* self) = v;
  //}  
  void Assign(PyObject* param) {
    /* note that these error does not raise error in python
       type check is actually done in wrapper layer */
    PyArrayObject *param0 = reinterpret_cast<PyArrayObject *>(param);
      
    if (!PyArray_Check(param0)){
       PyErr_SetString(PyExc_ValueError, "Input data must be ndarray");
       return;
    }
    int typ = PyArray_TYPE(param0);
    if (typ != NPY_DOUBLE){
        PyErr_SetString(PyExc_ValueError, "Input data must be float64");
	return;
    }
    int ndim = PyArray_NDIM(param0);
    if (ndim != 1){
      PyErr_SetString(PyExc_ValueError, "Input data NDIM must be one");
      return ;
    }
    npy_intp *shape = PyArray_DIMS(param0);    
    int len = self->Size();
    if (shape[0] != len){    
      PyErr_SetString(PyExc_ValueError, "input data length does not match");
      return ;
    }    
    (mfem::Vector &)(* self) = (double *) PyArray_DATA(param0);
  }
 
/*  this will be turn on in mfem-3.3.3
ParGridFunction(ParMesh *pmesh, const char *gf_file){
    mfem::ParGridFunction *pgf;
    mfem::ParFiniteElementSpace *pfes;    
    std::ifstream gfstream(gf_file);
    if (!gfstream)
    {
    std::cerr << "\nCan not open ParGridFunction: " << gf_file << '\n' << std::endl;
    return NULL;
    }
    pgf = new mfem::GridFunction(pmesh, gfstream);
    return pgf;
    }
*/ 
};  //end of extend
}  // end of namespace  

/*
   virtual void Save(std::ostream &out) const;
   void SaveAsOne(std::ostream &out = mfem::out);
*/
OSTREAM_ADD_DEFAULT_FILE(ParGridFunction, Save)
OSTREAM_ADD_DEFAULT_FILE(ParGridFunction, SaveAsOne)
