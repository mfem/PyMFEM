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

    dir = os.path.dirname(os.path.abspath(__file__))    
    meshfile = os.path.join(dir, "../data/amr-quad.mesh")
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()
    order = 2
    
    fec = mfem.DG_FECollection(order, dim)
    #fec = mfem.H1_FECollection(order, dim)    
    fes = mfem.FiniteElementSpace(mesh, fec, 1)
    print('Number of finite element unknowns: '+
          str(fes.GetTrueVSize()))
    print('Size of DoFVinite element unknowns: '+
          str(fes.GetNDofs()))

    points = np.zeros((fes.GetTrueVSize(), sdim))

    R = fes.GetConformingRestriction()
    if R is not None:    
        VDof2TDof = np.zeros(fes.GetNDofs(), dtype=int)
        for i, j in enumerate(R.GetJArray()):
            VDof2TDof[j] = i
        TDof2Vdof = R.GetJArray().copy()
    else:
        VDof2TDof = None
        TDof2VDof = None
        
    for j in range(fes.GetNE()):
        el = fes.GetFE(j)
        tr = fes.GetElementTransformation(j)
        vdofs = fes.GetElementVDofs(j)
        
        tdofs= vdofs if VDof2TDof is None else [VDof2TDof[k] for k in vdofs]
        
        ir = el.GetNodes()
        for k, tdof in enumerate(tdofs):
            points[tdof] = tr.Transform(ir.IntPoint(k))

    # search for tdof = 10
    tdof = 10
    vdof = tdof if TDof2VDof is None else TDof2Vdof[tdof]
    fes.BuildDofToArrays()
    j = fes.GetElementForDof(vdof)
    el = fes.GetFE(j)
    tr = fes.GetElementTransformation(j)
    vdofs = fes.GetElementVDofs(j)
    ir = el.GetNodes()
    for k, vv in enumerate(vdofs):
        if vdof == vv:
           pp = tr.Transform(ir.IntPoint(k))
           break
       
    #print(points)
    #print(pp)
if __name__=='__main__':
    run_test()
