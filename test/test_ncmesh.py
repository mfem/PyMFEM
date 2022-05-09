import sys
import os

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

def test_course_to_file_map1():
    dir = os.path.dirname(os.path.abspath(__file__))        
    mesh_file = os.path.join(dir, "../data/inline-quad.mesh")


    mesh = mfem.Mesh(mesh_file, 1, 1)

    refinements = mfem.RefinementArray()
    refinements.Append(mfem.Refinement(0, 0b11))
    refinements.Append(mfem.Refinement(1, 0b10))
    refinements.Append(mfem.Refinement(2, 0b01))

    mesh.GeneralRefinement(refinements)

    cft = mesh.GetRefinementTransforms()
    coarse_to_fine = mfem.Table()

    cft.MakeCoarseToFineTable(coarse_to_fine)
    print("coarse_to_fine element number mapping:")
    coarse_to_fine.Print()

    '''
    TODO: this part does not work due to API change.
    coarse_to_ref_type = mfem.intArray()
    ref_type_to_matrix = mfem.Table()
    ref_type_to_geom = mfem.GeometryTypeArray()

    cft.GetCoarseToFineMap(mesh,
                           coarse_to_fine,
                           coarse_to_ref_type,
                           ref_type_to_matrix,
                           ref_type_to_geom)

    print("ref_type_to_geom mapping:")

    for i in range(ref_type_to_geom.Size()):
        g = ref_type_to_geom[i]
        print("ref_type: " + str(i) + "=> geom: " + str(g))
   '''

def test_course_to_file_map2():
    dir = os.path.dirname(os.path.abspath(__file__))        
    mesh_file = os.path.join(dir, "../data/inline-quad.mesh")
    
    mesh = mfem.Mesh(mesh_file, 1, 1)

    els = mfem.intArray([0])
    mesh.GeneralRefinement(els)

    zero = mfem.Vector(mesh.GetNE())
    zero.Assign(0.0)

    print(zero.Size())
    mesh.DerefineByError(zero, 1.0)
    print(mesh.GetNE())
    dtrans = mesh.ncmesh.GetDerefinementTransforms()
    coarse_to_fine = mfem.Table()
    
    dtrans.GetCoarseToFineMap(mesh, coarse_to_fine)
    print("table rows: " + str(coarse_to_fine.Size()))
    
    for pe in range(coarse_to_fine.Size()):
        row = mfem.intArray()
        coarse_to_fine.GetRow(pe, row)
        nchild = row.Size();
        txt = "row " + str(pe) + ":"
        for fe in range(nchild):
            child = row[fe]
            txt += " " + str(child)
        print(txt)

if __name__ == '__main__':
    device = mfem.Device('cpu')
    
    test_course_to_file_map1()
    test_course_to_file_map2()

