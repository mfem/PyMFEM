%module nonlininteg
%{
#include "fem/nonlininteg.hpp"
#include "fem/linearform.hpp"  
#include "pycoefficient.hpp"
#include "pyoperator.hpp"     
%}
/*
%init %{
import_array();
%}
*/
%import vector.i
%import operators.i
%import fespace.i
%import eltrans.i
%include "fem/nonlininteg.hpp"
