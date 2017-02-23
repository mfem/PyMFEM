%module socketstream
%{
#include <iostream>
#include "iostream_typemap.hpp"      
#include "mesh/mesh_headers.hpp"
#include "fem/gridfunc.hpp"  
#include "general/socketstream.hpp"
#include "numpy/arrayobject.h"  
%}

%init %{
import_array();
%}
//%rename(sopen)  open(const char hostname[], int port);
%import "mesh.i"
%import "gridfunc.i"
%import "ostream_typemap.i"

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
  mfem::socketstream& __lshift__(const char ostr[])
   { 
      *self << ostr;
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
  mfem::socketstream& flush()
   { 
     self->flush();
   } 
} 
