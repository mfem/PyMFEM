%module(package="mfem._ser", directors="1")  lininteg
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "mfem.hpp"    
#include "pyoperator.hpp"      
#include "../common/pycoefficient.hpp"
#include "numpy/arrayobject.h"        
%}

%init %{
import_array();
%}

%include "exception.i"
 //%include "fem/coefficient.hpp"
%import "fe.i"
%import "vector.i"
%import "eltrans.i"
%import "intrules.i"
%import "coefficient.i"
%import "../common/exception_director.i"
%import "fespace.i"

%include "../common/lininteg_ext.i"

%include "fem/lininteg.hpp"


