from __future__ import print_function
import os
import sys

def run_test():
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
    t = mesh.GetEdgeVertexTable()
    t.Print()
    
if __name__=='__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':   
        import mfem.par as mfem
    else:
        import mfem.ser as mfem
        
    run_test()
