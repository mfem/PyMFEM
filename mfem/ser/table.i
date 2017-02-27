%module table
%{
#include "iostream_typemap.hpp"      
#include "general/table.hpp"
#include <iostream>
%}

%import "array.i"
%import "ostream_typemap.i"

%include "general/table.hpp"
