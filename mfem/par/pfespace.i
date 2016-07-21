%module pfespace
%{
#include <mpi.h>
#include "fem/pfespace.hpp"
#include "numpy/arrayobject.h"
#define MFEM_USE_MPI  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%inline %{
/*--------------------------------------------------------------------------
 * Big int stuff
 *--------------------------------------------------------------------------*/
#ifdef HYPRE_BIGINT
typedef long long int HYPRE_Int;
#define HYPRE_MPI_INT MPI_LONG_LONG_INT
#else 
typedef int HYPRE_Int;
#define HYPRE_MPI_INT MPI_INT
#endif
%}


%init %{
import_array();
%}

%import fespace.i
%import pmesh.i
%import hypre.i

%immutable face_nbr_glob_dof_map;
//DoF accesser

%feature("shadow") mfem::ParFiniteElementSpace::GetBdrElementDofs %{
def GetBdrElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _pfespace.ParFiniteElementSpace_GetBdrElementDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetElementDofs %{
def GetElementDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _pfespace.ParFiniteElementSpace_GetElementDofs(self, i, vdofs)
    return vdofs.ToList()
%}
%feature("shadow") mfem::ParFiniteElementSpace::GetFaceDofs %{
def GetFaceDofs(self, i):
    from  .array import intArray
    vdofs = intArray()
    _pfespace.ParFiniteElementSpace_GetFaceDofs(self, i, vdofs)
    return vdofs.ToList()
%}

#define MFEM_USE_MPI  
%include "fem/pfespace.hpp"
