//
// Copyright (c) 2020-2026, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") common_functions

%{
#include <iostream>
#include <fstream>
#include "mfem.hpp"
#include "../common/io_stream.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
  %}
%init %{
import_array1(-1);
%}

%include "exception.i"
%include "../common/exception.i"

%import "globals.i"
%import "array.i"
%import "globals.i"
%ignore "";
%rename("%s") mfem;
%rename("%s") mfem::Add;
%rename(Add) mfem::add;
%rename(Subtract) mfem::subtract;
%rename("%s") mfem::Transpose;
%rename("%s") mfem::Mult;
%rename("%s") mfem::RAP;
%rename("%s") mfem::InnerProduct;
%rename("%s") mfem::OuterProduct;
%rename("%s") mfem::TransposeMult;
%rename("%s") mfem::Mult_AtDA;
%newobject mfem::Add;
%newobject mfem::RAP;
%newobject mfem::Mult(const SparseMatrix&, DenseMatrix&);
%newobject mfem::Mult(const BlockMatrix&, const BlockMatrix&);
#ifdef MFEM_NUMBA_SWIG_BRIDGE_SER
%include "../common/bridges/common_functions_nsb_1.i"
#endif
%include "general/table.hpp"
%include "../common/vector_function_declarations.i"
%include "linalg/vector.hpp"
%include "linalg/densemat.hpp"
%include "linalg/sparsemat.hpp"
%include "linalg/blockmatrix.hpp"

#ifdef MFEM_NUMBA_SWIG_BRIDGE_SER
%include "../common/bridges/common_functions_nsb_2.i"
#endif
