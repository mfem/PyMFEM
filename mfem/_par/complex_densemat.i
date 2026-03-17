%module(package="mfem._par") complex_densemat
%{
#include "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/pyoperator.hpp"
#include "../common/pysolvers.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
#include "../common/pybilininteg.hpp"
#include "../common/pynonlininteg.hpp"
#include "../common/io_stream.hpp"
%}

%init %{
import_array();
%}

%include "exception.i"
%import "element.i"
%import "../common/exception.i"

%import "array.i"
%import "vector.i"
%import "densemat.i"
%import "complex_operator.i"
%import "../common/exception.i"
%import "../common/io_stream_typemap.i"

%ignore mfem::ComplexLUFactors::Mult(int m, int n, std::complex<real_t> *X) const;
%include "linalg/complex_densemat.hpp"
