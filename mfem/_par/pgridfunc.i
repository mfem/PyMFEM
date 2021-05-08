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
#include "../common/io_stream.hpp"      
#include  "config/config.hpp"
#include "mesh/pmesh.hpp"  
#include "fem/pgridfunc.hpp"
#include "fem/pfespace.hpp"  
#include "fem/linearform.hpp"  
#include "pycoefficient.hpp"  
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
   

%include "fem/pgridfunc.hpp"

namespace mfem{
%extend ParGridFunction{
ParGridFunction(mfem::ParFiniteElementSpace *fes, const mfem::Vector &v, int offset){
   mfem::ParGridFunction *gf;   
   gf = new mfem::ParGridFunction(fes, v.GetData() + offset);
   return gf;
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
