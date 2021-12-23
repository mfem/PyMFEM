%module(package="mfem._par", directors="0")  strumpack
%feature("autodoc", "1");

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
#include "../common/io_stream.hpp"         
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"
#include "pyoperator.hpp"

%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"
#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif


%include "exception.i"
%import "../common/numpy_int_typemap.i"

%include "../common/exception.i"
%import "operators.i"
%import "hypre.i"

// convert list input
 /*
%typemap(in) (int argc, char *argv[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = (char **) malloc(($1+1)*sizeof(char *));
  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
  $2[i] = 0;
}

%typemap(freearg) (int argc, char *argv[]) {
   if ($2) free($2);
}

%typecheck(SWIG_TYPECHECK_STRING_ARRAY) (int argc, char *argv[]) {
  $1 = PyList_Check($input) ? 1 : 0;
}
 */
%inline %{
char **argv_obj(PyObject* input){
  int i;
  if (!PyList_Check(input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  int num = PyList_Size(input);
  char **out = (char **) malloc((num+1)*sizeof(char *));
  for (i = 0; i < num; i++) {
    PyObject *s = PyList_GetItem(input,i);
    if (!PyString_Check(s)) {
        free(out);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    out[i] = PyString_AsString(s);
  }
  out[i] = 0;
  return out;
 };
 %}

%newobject argv_obj;
//%include "carrays.i"

//%array_class(char *, ptcharArray);
//%pointer_functions(char, charp);

%pythonprepend mfem::STRUMPACKSolver::STRUMPACKSolver %{
  attach_argv = False
  if isinstance(args[0], list):
      aa = [""]+args[0]  
      num = len(aa)
      ptr = argv_obj(aa)
      args = (num, ptr, args[1])
      attach_argv = True

%}

%pythonappend mfem::STRUMPACKSolver::STRUMPACKSolver %{
  if attach_argv:
      self._argv = ptr
%}  

%include "SPOptions.hpp"
%include "linalg/strumpack.hpp"

