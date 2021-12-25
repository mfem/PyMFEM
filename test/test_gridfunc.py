from __future__ import print_function
import os
from os.path import expanduser, join
import sys
import numpy as np
import io
from mfem import path as mfem_path

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
    #meshfile = expanduser(join(mfem_path, 'data', 'semi_circle.mesh'))
    dir = os.path.dirname(os.path.abspath(__file__))
    meshfile = os.path.join(dir, "../data/amr-quad.mesh")
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()    
    fec = mfem.H1_FECollection(1, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec, 1)
    print('Number of finite element unknowns: '+
          str(fespace.GetTrueVSize()))
    
    c = mfem.ConstantCoefficient(1.0)
    gf = mfem.GridFunction(fespace)
    gf.ProjectCoefficient(c)

    print("write mesh to STDOUT")
    mesh.Print(mfem.STDOUT)
    print("creat VTK file to file")
    mesh.PrintVTK('mesh.vtk', 1)
    print("creat VTK to STDOUT")
    mesh.PrintVTK(mfem.STDOUT, 1)
    print("save GridFunction to file")
    gf.Save('out_test_gridfunc.gf')
    gf.SaveVTK(mfem.wFILE('out_test_gridfunc1.vtk'), 'data', 1)
    print("save GridFunction to file in VTK format")    
    gf.SaveVTK('out_test_gridfunc2.vtk', 'data', 1)
    print("Gridfunction to STDOUT")
    gf.Save(mfem.STDOUT)
    
    o = io.StringIO()
    count = gf.Save(o)
    count2 = gf.SaveVTK(o, 'data', 1)
    print("length of data ", count, count2)
    print('result: ', o.getvalue())    
    
if __name__=='__main__':
    run_test()
