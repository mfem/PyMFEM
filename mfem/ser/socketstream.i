%module socketstream
%{
#include <iostream>  
#include "general/socketstream.hpp"
%}
%rename(sopen)  open(const char hostname[], int port);
%include "general/socketstream.hpp"
