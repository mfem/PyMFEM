%module pgridfunc
%{
#include <mpi.h>
#include "iostream_typemap.hpp"      
#include  "config/config.hpp"
#include "fem/pgridfunc.hpp"
#include "fem/linearform.hpp"  
#include "pycoefficient.hpp"  
#include "numpy/arrayobject.h"
%}
%include  "config/_config.hpp" // include mfem MACRO

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}

%include "../common/cpointers.i"
%import pfespace.i
%import gridfunc.i
%import hypre.i
%import pmesh.i
%import linearform.i

%import "ostream_typemap.i"

%pointer_class(int, intp);

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

%rename(Assign) mfem::ParGridFunction::operator=;

%newobject mfem::ParGridFunction::ParallelAssemble;
%newobject mfem::ParGridFunction::ParallelAverage;
%newobject mfem::ParGridFunction::ParallelProject;
%newobject mfem::ParGridFunction::GetTrueDofs;

%include "fem/pgridfunc.hpp"

namespace mfem{
%extend ParGridFunction{
ParGridFunction(mfem::ParFiniteElementSpace *fes, const mfem::Vector &v, int offset){
   mfem::ParGridFunction *gf;   
   gf = new mfem::ParGridFunction(fes, v.GetData() + offset);
   return gf;
 }
   }  //end of extend
}  // end of namespace  
