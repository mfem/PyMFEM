%module pmesh
%{
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>  
#include <mpi.h>
#include "iostream_typemap.hpp"     
#include "config/config.hpp"
#include "mesh/pmesh.hpp"
#include "fem/linearform.hpp"
#include "general/communication.hpp"  
#include "numpy/arrayobject.h"
%}

%include  "config/_config.hpp" // include mfem MACRO
%init %{
import_array();
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%include "../common/cpointers.i"
%import mesh.i
%import pncmesh.i
 //
%import communication.i
%import "ostream_typemap.i"

%immutable face_nbr_elements;
%immutable face_nbr_vertices;
%immutable gtopo;

%pointer_class(int, intp);

%feature("shadow") mfem::ParMesh::GroupFace %{
def GroupFace(self, group, i, *args):
    if len(args) == 0:
        args = (face = intp()
        o = intp()
        _pmesh.Mesh_GroupFace(self, group, i, edge, o)      
        return edge.value(), o.value()
    else:
        _pmesh.Mesh_GroupFace(self, group, i, *args)            
%}
%feature("shadow") mfem::ParMesh::GroupEdge ${
def GroupEdge(self, group, i, *args):
    if len(args) == 0:
        edge = intp()
        o = intp()  
        _pmesh.Mesh_GroupEdge(self, group, i, edge, o)
        return edge.value(), o.value()
    else:
        _pmesh.Mesh_GroupEdge(self, group, i, *args)      
%}


%include "mesh/pmesh.hpp"

namespace mfem{
%extend ParMesh{
ParMesh(MPI_Comm comm, const char *mesh_file){
    mfem::ParMesh *mesh;
    std::ifstream imesh(mesh_file);
    if (!imesh)
    {
    std::cerr << "\nCan not open mesh file: " << mesh_file << '\n' << std::endl;
    return NULL;
    }
    mesh = new mfem::ParMesh(comm, imesh);
    return mesh;
    }
void ParPrintToFile(const char *mesh_file, const int precision) const
    {
    std::ofstream mesh_ofs(mesh_file);	
    mesh_ofs.precision(precision);
    self->ParPrint(mesh_ofs);	
    }
};   
}

