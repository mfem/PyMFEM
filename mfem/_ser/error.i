%module(package="mfem._ser") error
%{
#include <exception>
#include <iostream>
#include <stdio.h>
#include "general/error.hpp"
%}
%include "exception.i"
%include "../common/exception.i"

%include "general/error.hpp"
