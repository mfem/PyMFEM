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
    myid = MPI.COMM_WORLD.rank
    smyid = '.'+'{:0>6d}'.format(myid)
else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0
    smyid = ''


def check(a, b, msg):
    assert len(a) == len(b), msg
    assert np.sum(np.abs(a - b)) == 0, msg


def check_mesh(m1, m2, msg):
    assert m1.GetNE() == m2.GetNE(), msg
    assert m1.GetNBE() == m2.GetNBE(), msg
    assert m1.GetNV() == m2.GetNV(), msg


def run_test():
    #meshfile = expanduser(join(mfem_path, 'data', 'semi_circle.mesh'))
    dir = os.path.dirname(os.path.abspath(__file__))    
    meshfile = os.path.join(dir, "../data/amr-quad.mesh")
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()
    fec = mfem.H1_FECollection(1, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec, 1)
    print('Number of finite element unknowns: ' +
          str(fespace.GetTrueVSize()))

    c = mfem.ConstantCoefficient(1.0)
    gf = mfem.GridFunction(fespace)
    gf.ProjectCoefficient(c)

    odata = gf.GetDataArray().copy()

    gf.Save("out_test_gz"+smyid+".gf")
    gf2 = mfem.GridFunction(mesh, "out_test_gz"+smyid+".gf")
    odata2 = gf2.GetDataArray().copy()
    check(odata, odata2, "text file does not agree with original")

    gf.Save("out_test_gz"+smyid+".gz")
    gf2.Assign(0.0)
    gf2 = mfem.GridFunction(mesh, "out_test_gz"+smyid+".gz")
    odata2 = gf2.GetDataArray().copy()
    check(odata, odata2, ".gz file does not agree with original")

    gf.Print("out_test_gz"+smyid+".dat")
    gf2.Assign(0.0)
    gf2.Load("out_test_gz"+smyid+".dat", gf.Size())
    odata2 = gf2.GetDataArray().copy()
    check(odata, odata2, ".dat file does not agree with original")

    gf.Print("out_test_gz.dat"+smyid+".gz")
    gf2.Assign(0.0)
    gf2.Load("out_test_gz.dat"+smyid+".gz", gf.Size())
    odata2 = gf2.GetDataArray().copy()
    check(odata, odata2, ".dat file does not agree with original (gz)")

    import gzip
    import io
    gf.Print("out_test_gz.dat2"+smyid+".gz")
    with gzip.open("out_test_gz.dat2"+smyid+".gz", 'rt') as f:
        sio = io.StringIO(f.read())
    gf3 = mfem.GridFunction(fespace)
    gf3.Load(sio, gf.Size())
    odata3 = gf3.GetDataArray().copy()
    check(odata, odata3, ".dat file does not agree with original(gz-io)")

    c = mfem.ConstantCoefficient(2.0)
    gf.ProjectCoefficient(c)
    odata = gf.GetDataArray().copy()

    o = io.StringIO()
    gf.Print(o)
    gf2.Load(o, gf.Size())
    odata2 = gf2.GetDataArray().copy()
    check(odata, odata2, "StringIO does not agree with original")

    print("GridFunction .gf, .gz .dat and StringIO agree with original")

    mesh2 = mfem.Mesh()
    mesh.Print("out_test_gz.mesh")
    mesh2.Load("out_test_gz.mesh")
    check_mesh(mesh, mesh2, ".mesh does not agree with original")

    mesh2 = mfem.Mesh()
    mesh.Print("out_test_gz.mesh"+smyid+".gz")
    mesh2.Load("out_test_gz.mesh"+smyid+".gz")

    check_mesh(mesh, mesh2, ".mesh.gz does not agree with original")

    mesh3 = mfem.Mesh()    
    mesh.PrintGZ("out_test_gz3"+smyid+".mesh")
    mesh3.Load("out_test_gz3"+smyid+".mesh")


    check_mesh(mesh, mesh3, ".mesh (w/o .gz exntension) does not agree with original")
    

    o = io.StringIO()        
    mesh2 = mfem.Mesh()    

    mesh.Print(o)
    mesh2.Load(o)
    check_mesh(mesh, mesh2, ".mesh.gz does not agree with original")

    print("Mesh .mesh, .mesh.gz and StringIO agree with original")

    print("PASSED")


if __name__ == '__main__':
    run_test()
