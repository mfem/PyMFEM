%module(package="mfem._ser") symmat
%{
#include  "mfem.hpp"
#include "linalg/symmat.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"    
#include "pyoperator.hpp"    
%}

%init %{
import_array();
%}
%include "exception.i"

%import "globals.i"
%import "matrix.i"
%include "../common/exception.i"


%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)

%include "linalg/symmat.hpp"

OSTREAM_ADD_DEFAULT_FILE(DenseSymmetricMatrix, Print)
