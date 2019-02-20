%module(directors="0")  pumi
%{
#include "mesh/mesh_headers.hpp"
#include "pumi.hpp"
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

%}

%init %{
import_array();
%}

%include "../common/mfem_config.i"
%include "../common/cpointers.i"
%include "exception.i"
%import "ostream_typemap.i"
%import "../common/numpy_int_typemap.i"
%import "mesh.i"
%import "pmesh.i"

%include "../common/exception.i"
%include "mesh/pumi.hpp"


