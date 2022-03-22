from __future__ import print_function
import os
import sys
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == '-p':   
    import mfem.par as mfem
    use_parallel = True
    from mfem.common.mpi_debug import nicePrint as print
    from mpi4py import MPI
    myid  = MPI.COMM_WORLD.rank
    
else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0

def InArray(T, sz, i):
    hoge = mfem.intArray((T, sz))
    ret = (i in hoge.ToList())
    return ret
    
def IndicesAreConnected(e2e, i, j):
    return (InArray(e2e.GetRow(i), e2e.RowSize(i), j) and
            InArray(e2e.GetRow(j), e2e.RowSize(j), i))

def run_test():
    n = 3

    ### 1D Periodic mesh
    orig_mesh = mfem.Mesh_MakeCartesian1D(n)
    translations = (mfem.Vector([1.0]),)
    #vv =mfem.vector_Vector(translations)
    mapping = orig_mesh.CreatePeriodicVertexMapping(translations)

    mesh = mfem.Mesh_MakePeriodic(orig_mesh, mapping)

    assert mesh.GetNV() == n, "GetNV is not correct"

    e2e = mesh.ElementToElementTable()
    assert IndicesAreConnected(e2e, 0, 2), "wrong connection"
    assert IndicesAreConnected(e2e, 0, 1), "wrong connection"
    assert IndicesAreConnected(e2e, 1, 2), "wrong connection"

    ### 2D Periodic mesh
    els = (mfem.Element.TRIANGLE, mfem.Element.QUADRILATERAL)

    for el in els:
        orig_mesh = mfem.Mesh_MakeCartesian2D(n, n, el, False, 1, 1, False)
        translations = (mfem.Vector([1., 0.]),
                        mfem.Vector([0., 1.]),)    
        mapping = orig_mesh.CreatePeriodicVertexMapping(translations)
        mesh = mfem.Mesh_MakePeriodic(orig_mesh, mapping)

        assert mesh.GetNV() == (n-1)**2 + 2*(n-1) + 1, "GetNV is not correct"

        if el == mfem.Element.QUADRILATERAL:
            e2e = mesh.ElementToElementTable()
            for i in range(n):
                assert IndicesAreConnected(e2e, i, i + n*(n-1)), "wrong connection"
                assert IndicesAreConnected(e2e, i*n, n-1 + i*n), "wrong connection"

    ### 3D Periodic mesh
    els = (mfem.Element.TETRAHEDRON,
           mfem.Element.HEXAHEDRON,
           mfem.Element.WEDGE)
    for el in els:
        orig_mesh = mfem.Mesh_MakeCartesian3D(n, n, n, el, 1, 1, 1, False)
        translations = (mfem.Vector([1., 0., 0.]),
                        mfem.Vector([0., 1., 0.]),
                        mfem.Vector([0., 0., 1.]),)
        mapping = orig_mesh.CreatePeriodicVertexMapping(translations)
        mesh = mfem.Mesh_MakePeriodic(orig_mesh, mapping)

        assert mesh.GetNV() == (n-1)**3 + 3*(n-1)**2 + 3*(n-1) +1, "GetNV is not correct"        
        if el == mfem.Element.HEXAHEDRON:
            e2e = mesh.ElementToElementTable()
            n2 = n*n
            for j in range(n):
                for i in range(n):       
                    assert IndicesAreConnected(e2e,
                                               i + j*n,
                                               i + j*n + n2*(n-1)), "wrong connection"
                    assert IndicesAreConnected(e2e,
                                               i + j*n2,
                                               i + j*n2 + n*(n-1)), "wrong connection"
                    assert IndicesAreConnected(e2e,
                                               i*n + j*n2,
                                               i*n + j*n2 + n-1), "wrong connection"



if __name__=='__main__':
    run_test()


