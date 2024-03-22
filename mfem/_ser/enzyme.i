%module(package="mfem._ser") enzyme
%{
#include "mfem.hpp"
%}

%include "../common/existing_mfem_headers.i"
#ifdef FILE_EXISTS_GENERAL_ENZYME

%import "exception.i"

%include "general/enzyme.hpp"

#endif
