//
// Copyright (c) 2020-2026, Princeton Plasma Physics Laboratory, All rights reserved.
//
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
