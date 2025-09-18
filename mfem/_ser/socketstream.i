//
// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.
//
%module(package="mfem._ser") socketstream
%feature("autodoc", "1");

%{
#include <iostream>
#include  "mfem.hpp"
#include "numpy/arrayobject.h"
#include "../common/io_stream.hpp"
#include "../common/pyoperator.hpp"
#include "../common/pycoefficient.hpp"
#include "../common/pyintrules.hpp"
%}

%init %{
import_array();
%}

//%rename(sopen)  open(const char hostname[], int port);
%include "exception.i"
%import "mesh.i"
%import "gridfunc.i"
%import "../common/exception.i"

%include "general/socketstream.hpp"

%extend mfem::socketstream{
  int precision(const int p)
   {
     return self->precision(p);
   }
  int precision()
   {
     return self->precision();
   }
  void send_solution(const mfem::Mesh &mesh,
 				    const mfem::GridFunction &gf)
   {
     *self << "solution\n" << mesh << gf << std::flush;
   }
  void send_text(const char ostr[])
   {
      *self << ostr << std::endl;
   }
  void flush()
   {
     *self << std::flush;
   }
  bool good()
  {
    return self->std::iostream::good();
  }

  mfem::socketstream& __lshift__(const char ostr[])
   {
      *self << ostr;
      return *self;
   }
  mfem::socketstream& __lshift__(const int x)
   {
      *self << x;
      return *self;
   }
  mfem::socketstream& __lshift__(const mfem::Mesh &mesh)
   {
      *self << mesh;
      return *self;
   }
  mfem::socketstream& __lshift__(const mfem::GridFunction &gf)
   {
      *self << gf;
      return *self;
   }
  mfem::socketstream& endline()
   {
     *self << std::endl;
     return *self;
   }


}
