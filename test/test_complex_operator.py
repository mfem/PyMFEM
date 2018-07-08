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

    if use_parallel:
        mesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)        
        fes = mfem.ParFiniteElementSpace(mesh, fec)        
        a1 = mfem.ParBilinearForm(fes)
        a2 = mfem.ParBilinearForm(fes)    
    else:
        fes = mfem.FiniteElementSpace(mesh, fec)        
        a1 = mfem.BilinearForm(fes)
        a2 = mfem.BilinearForm(fes)    
    one = mfem.ConstantCoefficient(1.0)    
    a1.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
    a1.Assemble()
    a1.Finalize();

    a2.AddDomainIntegrator(mfem.DiffusionIntegrator(one))
    a2.Assemble()
    a2.Finalize()
    
    if use_parallel:
        M1 = a1.ParallelAssemble()
        M2 = a2.ParallelAssemble()
        M1.Print('M1')
        width = fes.GetTrueVSize()
        #X = mfem.HypreParVector(fes)
        #Y = mfem.HypreParVector(fes)
        #X.SetSize(fes.TrueVSize())
        #Y.SetSize(fes.TrueVSize())
        #from mfem.common.parcsr_extra import ToScipyCoo
        #MM1 = ToScipyCoo(M1)
        #print(MM1.toarray())
        #print(MM1.dot(np.ones(6)))        
    else:
        M1 = a1.SpMat()        
        M2 = a2.SpMat()
        M1.Print('M1')
        width = fes.GetVSize()
        
        #X = mfem.Vector()
        #Y = mfem.Vector()
        #X.SetSize(M1.Width())
        #Y.SetSize(M1.Height())
        #from mfem.common.sparse_utils import sparsemat_to_scipycsr
        #MM1 = sparsemat_to_scipycsr(M1, np.float)
        #print(MM1.toarray())        
        #print(MM1.dot(np.ones(6)))
    #X.Assign(0.0)
    #X[0] = 1.0
    #M1.Mult(X, Y)
    #print(Y.GetDataArray())
    
    Mc = mfem.ComplexOperator(M1, M2, hermitan=True)
    offsets = mfem.intArray([0, width, width])
    offsets.PartialSum()

    x = mfem.BlockVector(offsets)
    y = mfem.BlockVector(offsets)

    x.GetBlock(0).Assign(0)
    if myid == 0:
        x.GetBlock(0)[0] = 1.0
    x.GetBlock(1).Assign(0)
    if myid == 0:
        x.GetBlock(1)[0] = 1.0
    
    Mc.Mult(x, y)
    print("x", x.GetDataArray())    
    print("y", y.GetDataArray())
    
    if myid == 0:
        x.GetBlock(1)[0] = -1.0

    x.Print()    
    Mc.Mult(x, y)
    print("x", x.GetDataArray())    
    print("y", y.GetDataArray())
    
if __name__=='__main__':
    run_test()
