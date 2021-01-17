from __future__ import print_function
import os
import sys

def run_test1(mfem):
    print("Test mesh module")
    Nvert = 6; Nelem = 8; Nbelem=2
    
    mesh = mfem.Mesh(2, Nvert, Nelem, 2, 3)
    tri_v = [[1.,  0.,  0.], [0.,  1.,  0.], [-1.,  0.,  0.],
             [0., -1.,  0.], [0.,  0.,  1.], [ 0.,  0., -1.]]
    tri_e = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4],
             [1, 0, 5], [2, 1, 5], [3, 2, 5], [0, 3, 5]]
    tri_l = [[1,4], [1,2]]
    
    for j in range(Nvert):
        mesh.AddVertex(tri_v[j])
    for j in range(Nelem):
        mesh.AddTriangle(tri_e[j], 1)
    for j in range(Nbelem):
        mesh.AddBdrSegment(tri_l[j], 1)
        
    mesh.FinalizeTriMesh(1,1, True)

    print(mesh.GetEdgeVertices(1))
    print(mesh.GetFaceElements(1))
    
    assert mesh.GetBdrElementAdjacentElement(1) == (1, 0), "GetBdrElementadjacentElement failed"

    i = 0
    Tr = mesh.GetElementTransformation(i)
    print(Tr.Weight())
    Tr = mesh.GetBdrElementTransformation(i)
    print(Tr.Weight())    
    Tr = mesh.GetFaceTransformation(i)
    print(Tr.Weight())    
    Tr = mesh.GetEdgeTransformation(i)
    print(Tr.Weight())

def run_test2(mfem):
    mesh = mfem.Mesh(6, 6, 6, "TETRAHEDRON")
    #r = mesh.CartesianPartitioning([2, 2, 2])
    #print(r)

    d = mfem.intArray([2, 2, 2])
    dd = d.GetData()

    r = mesh.CartesianPartitioning(dd)
    result = mfem.intArray()
    result.MakeRef(r, mesh.GetNE())
    result.MakeDataOwner()
    print(result.ToList())

    '''    

    mesh.Print("sample.mesh")
    r = mesh.CartesianPartitioning([2, 2, 2])
    print(r)
    '''
if __name__=='__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':   
        import mfem.par as mfem
    else:
        import mfem.ser as mfem
        
    run_test1(mfem)
    run_test2(mfem)    
