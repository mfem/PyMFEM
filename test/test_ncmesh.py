import sys

if len(sys.argv) > 1 and sys.argv[1] == '-p':
    import mfem.par as mfem
    use_parallel = True

    from mfem.common.mpi_debug import nicePrint as print

    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0

def test_course_to_file_map():
    mesh_file = "../data/inline-quad.mesh"
    device = mfem.Device('cpu')

    mesh = mfem.Mesh(mesh_file, 1, 1)

    refinements = mfem.RefinementArray()
    refinements.Append(mfem.Refinement(0, 0b11))
    refinements.Append(mfem.Refinement(1, 0b10))
    refinements.Append(mfem.Refinement(2, 0b01))

    mesh.GeneralRefinement(refinements)

    cft = mesh.GetRefinementTransforms()

    coarse_to_fine = mfem.Table()
    coarse_to_ref_type = mfem.intArray()
    ref_type_to_matrix = mfem.Table()
    ref_type_to_geom = mfem.GeometryTypeArray()

    cft.GetCoarseToFineMap(mesh,
                           coarse_to_fine,
                           coarse_to_ref_type,
                           ref_type_to_matrix,
                           ref_type_to_geom)
    print("coarse_to_fine element number mapping:")
    coarse_to_fine.Print()

    print("ref_type_to_geom mapping:")

    for i in range(ref_type_to_geom.Size()):
        g = ref_type_to_geom[i]
        print("ref_type: " + str(i) + "=> geom: " + str(g))

if __name__ == '__main__':
    test_course_to_file_map()

'''
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   const char *mesh_file = "../data/inline-quad.mesh";

   const char *device_config = "cpu";
   Device device(device_config);

   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   Array<Refinement> refinements;
   refinements.Append(Refinement(0,0b11));
   refinements.Append(Refinement(1,0b10));
   refinements.Append(Refinement(2,0b01));
 
   mesh.GeneralRefinement(refinements);

   mfem::CoarseFineTransformations const & coarse_fine_transform
      = mesh.GetRefinementTransforms();
   mfem::Table coarse_to_fine;
   mfem::Array<int> coarse_to_ref_type;
   mfem::Table ref_type_to_matrix;
   mfem::Array<mfem::Geometry::Type> ref_type_to_geom;
   coarse_fine_transform.GetCoarseToFineMap(mesh, coarse_to_fine,
                                            coarse_to_ref_type,
                                            ref_type_to_matrix,
                                            ref_type_to_geom);

   printf("coarse_to_fine element number mapping:\n");
   coarse_to_fine.Print();

   printf("ref_type_to_geom mapping:\n");
   for (int i = 0; i < ref_type_to_geom.Size(); i++) {
      Geometry::Type g = ref_type_to_geom[i];
      printf("ref_type: %d => geom: %d\n",i,g);
   }
}
'''
