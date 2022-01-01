%module(package="mfem._par") pmesh

%feature("autodoc", "1");

%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>  
#include <mpi.h>
#include "mfem.hpp"
#include "pyoperator.hpp"
#include "../common/io_stream.hpp"     
#include "numpy/arrayobject.h"
#include "../common/pycoefficient.hpp"    
%}

%include "../common/mfem_config.i"

%init %{
import_array();
%}

#ifdef MFEM_USE_MPI
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%include "exception.i"
%include "std_string.i"

 //%include "../common/cpointers.i"
 //%import "cpointers.i"
%import "mesh.i"
%import "pncmesh.i"
%import "hypre.i"
%import "communication.i"
%import "../common/exception.i"

%import "../common/io_stream_typemap.i"
OSTREAM_TYPEMAP(std::ostream&)
ISTREAM_TYPEMAP(std::istream&)

%immutable face_nbr_elements;
%immutable face_nbr_vertices;
%immutable gtopo;

%pointer_class(int, intp);

%feature("shadow") mfem::ParMesh::GroupFace %{
def GroupFace(self, group, i, *args):
    if len(args) == 0:
        from mfem.par import intp    
        face = intp()
        o = intp()
        $action(self, group, i, face, o)      
        return face.value(), o.value()
    else:
        return $action(self, group, i, *args)            
%}
	  
%feature("shadow") mfem::ParMesh::GroupEdge %{
def GroupEdge(self, group, i, *args):
    if len(args) == 0:
        from mfem.par import intp  
        edge = intp()
        o = intp()  
        $action(self, group, i, edge, o)
        return edge.value(), o.value()
    else:
        return $action(self, group, i, *args)      
%}


%include "mesh/pmesh.hpp"

namespace mfem{
%extend ParMesh{
     //     
     //ParMesh(MPI_Comm comm, const char *mesh_file){
     //mfem::ParMesh *mesh;
     //std::ifstream imesh(mesh_file);
     //if (!imesh)
     //{
     //std::cerr << "\nCan not open mesh file: " << mesh_file << '\n' << std::endl;
     //return NULL;
     //}
     //mesh = new mfem::ParMesh(comm, imesh);
     //return mesh;
     //}
void ParPrintToFile(const char *mesh_file, const int precision) const
    {
    std::ofstream mesh_ofs(mesh_file);	
    mesh_ofs.precision(precision);
    self->ParPrint(mesh_ofs);	
    }
};   
}

/*
  virtual void Print(std::ostream &out = mfem::out) const;
  virtual void PrintXG(std::ostream &out = mfem::out) const;
  void PrintAsOne(std::ostream &out = mfem::out);
  void PrintAsOneXG(std::ostream &out = mfem::out);
  virtual void PrintInfo(std::ostream &out = mfem::out);
  void ParPrint(std::ostream &out) const;
*/

#ifndef SWIGIMPORTED
OSTREAM_ADD_DEFAULT_FILE(ParMesh, Print)
OSTREAM_ADD_DEFAULT_FILE(ParMesh, PrintXG)
OSTREAM_ADD_DEFAULT_FILE(ParMesh, PrintAsOne)
OSTREAM_ADD_DEFAULT_FILE(ParMesh, PrintAsOneXG)
OSTREAM_ADD_DEFAULT_FILE(ParMesh, PrintInfo)
OSTREAM_ADD_DEFAULT_STDOUT_FILE(ParMesh, ParPrint)
#endif
