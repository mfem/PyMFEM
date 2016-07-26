%module nonlinearform
%{
#include "fem/nonlininteg.hpp"
#include "fem/nonlinearform.hpp"
#include "pyoperator.hpp"     
%}
/*
%init %{
import_array();
%}
*/
%import operators.i
%import fespace.i
%import nonlininteg.i
%include "fem/nonlinearform.hpp"
