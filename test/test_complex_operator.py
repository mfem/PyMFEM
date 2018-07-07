from __future__ import print_function
import os
import sys

if len(sys.argv) > 1 and sys.argv[1] == '-p':   
    import mfem.par as mfem
    use_parallel = True
else:
    import mfem.ser as mfem
    use_parallel = False        

def run_test():
    print("Test complex_operator module")
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
    dim = mesh.Dimension()
    order = 1
    fec = mfem.H1_FECollection(order, dim)
    fes = mfem.FiniteElementSpace(mesh, fec)
    
    a1 = mfem.BilinearForm(fes)
    a2 = mfem.BilinearForm(fes)    
    one = mfem.ConstantCoefficient(1.0)    
    a1.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
    a1.Assemble()
    a1.Finalize();
    M1 = a1.SpMat()
    a2.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
    a2.Assemble()
    a2.Finalize();
    M2 = a2.SpMat()
    Mc = mfem.ComplexOperator(M1, M2, hermitan=True)
    
if __name__=='__main__':
    run_test()
