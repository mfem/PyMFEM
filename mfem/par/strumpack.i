%module(directors="0")  strumpack
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
#include "iostream_typemap.hpp"         
#include "numpy/arrayobject.h"
#include "pycoefficient.hpp"
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


%include "../common/cpointers.i"
%include "exception.i"
%import "ostream_typemap.i"
%import "../common/numpy_int_typemap.i"

%include "../common/exception.i"
%import "operators.i"
%import "hypre.i"

%include "SPOptions.hpp"
%include "linalg/strumpack.hpp"


