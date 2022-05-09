from __future__ import print_function

import io
import gzip

import numpy as np

import os
from os.path import expanduser, join
import sys
import numpy as np
import io

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

def check(a, b, msg):
    assert len(a) == len(b), msg
    assert np.sum(np.abs(a-b)) == 0, msg

def check_mesh(m1, m2, msg):
    assert m1.GetNE() == m2.GetNE(), msg
    assert m1.GetNBE() == m2.GetNBE(), msg
    assert m1.GetNV() == m2.GetNV(), msg

def run_test():    
    mesh = mfem.Mesh(10, 10, "TRIANGLE")

    fec = mfem.H1_FECollection(1, mesh.Dimension())
    fespace = mfem.FiniteElementSpace(mesh, fec)

    gf = mfem.GridFunction(fespace)
    gf.Assign(1.0)
    odata = gf.GetDataArray().copy()

    visit_dc = mfem.VisItDataCollection("test_gf", mesh)
    visit_dc.RegisterField("gf", gf)
    visit_dc.Save()

    mesh2 = mfem.Mesh("test_gf_000000/mesh.000000")
    check_mesh(mesh, mesh2, ".mesh2 does not agree with original")

    gf2 = mfem.GridFunction(mesh2, "test_gf_000000/gf.000000")
    odata2 = gf2.GetDataArray().copy()
    check(odata, odata2, "odata2 file does not agree with original")

    visit_dc = mfem.VisItDataCollection("test_gf_gz", mesh)
    visit_dc.SetCompression(True)
    visit_dc.RegisterField("gf", gf)
    visit_dc.Save()

    with gzip.open("test_gf_gz_000000/mesh.000000", 'rt') as f:
       sio = io.StringIO(f.read())
    mesh3 = mfem.Mesh(sio)
    check_mesh(mesh, mesh3, ".mesh3 does not agree with original")

    with gzip.open("test_gf_gz_000000/gf.000000", 'rt') as f:
       sio = io.StringIO(f.read())

    # This is where the error is:
    gf3 = mfem.GridFunction(mesh3, sio)
    odata3 = gf3.GetDataArray()
    check(odata, odata3, "gf3 file does not agree with original")

    mesh4 = mfem.Mesh("test_gf_gz_000000/mesh.000000")
    check_mesh(mesh, mesh4, ".mesh4 does not agree with original")

    gf4 = mfem.GridFunction(mesh3, "test_gf_gz_000000/gf.000000")
    odata4 = gf4.GetDataArray()

    check(odata, odata4, "gf3 file does not agree with original")

    print("PASSED")    

if __name__=='__main__':
    run_test()
