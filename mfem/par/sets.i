%module sets
%{
#include <mpi.h>
#include "general/sets.hpp"
%}
/*
%init %{
import_array();
%}
*/
%import array.i
%import table.i
%include "general/sets.hpp"
