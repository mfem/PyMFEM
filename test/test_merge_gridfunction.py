#
# This script tests an approach to merge ParGridFunction to GridFunction
#
from __future__ import print_function
import os
from os.path import expanduser, join
import sys
import numpy as np
import io
from mfem import path as mfem_path

import mfem.par as mfem
use_parallel = True
from mpi4py import MPI

from mfem.common.mpi_debug import nicePrint

@mfem.jit.scalar()
def func(ptx):
    return ptx[0]

fec_type = mfem.H1_FECollection
fec_order = 1

myid = MPI.COMM_WORLD.rank
nprc = MPI.COMM_WORLD.size

def generate_data():
    ### generate serial mesh
    smesh = mfem.Mesh(3, 3, 3, "TETRAHEDRON")
    #smesh.ReorientTetMesh()
    
    dim = smesh.Dimension()
    sdim = smesh.SpaceDimension()
    
    fec = fec_type(fec_order, dim)    
    fes = mfem.FiniteElementSpace(smesh, fec, 1)
    
    ### make a reference GF, simply casting a coefficient 
    x_ser = mfem.GridFunction(fes)
    x_ser.ProjectCoefficient(func)


    ### Generate patitioning, we need this to merge solutiosn later
    parts = smesh.GeneratePartitioning(MPI.COMM_WORLD.size, 1)

    mesh = mfem.Mesh(3, 3, 3, "TETRAHEDRON")
    #mesh.ReorientTetMesh()
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh, parts)
    
    pfes = mfem.ParFiniteElementSpace(pmesh, fec, 1)    
    x_par = mfem.ParGridFunction(pfes)    
    x_par.ProjectCoefficient(func)

    ### Save data set to files
    smyid = '.{:0>6d}'.format(myid)        
    x_par.Save('merge_sol_data'+smyid, 16)
    pmesh.Print('merge_sol_mesh'+smyid)
    if myid == 0:
        smesh.Save("merge_sol.mesh")
        x_ser.Save("merge_sol_ref.gf")

    return parts

def merge_data(parts):
    ### load mesh and create FES, then collect VDofs for all elements
    smesh = mfem.Mesh("merge_sol.mesh")

    fec = fec_type(fec_order, smesh.Dimension())    
    fes = mfem.FiniteElementSpace(smesh, fec, 1)

    # convert parts (* int) to intArray so that we can iterate over it.
    idx = mfem.intArray((parts, smesh.GetNE()))

    x_merge = mfem.GridFunction(fes)
    x_merge.Assign(0)
        
    for myid2 in range(nprc):
        smyid2 = '.{:0>6d}'.format(myid2)
        m = mfem.Mesh('merge_sol_mesh'+smyid2)
        #m = mfem.Mesh('merge_sol_mesh'+smyid2, 0, 0, False) # refine flag should be 1
        #m.ReorientTetMesh()                        
        gf = mfem.GridFunction(m, 'merge_sol_data'+smyid2)
        fes2 = gf.FESpace()

        iel = 0
        for k, rank in enumerate(idx):
            if myid2 != rank:
                continue

            svdofs = fes.GetElementVDofs(k)    # this is VDofs from serial mesh
            vdofs2 = fes2.GetElementVDofs(iel) # this is VDofs from parallel mesh
            for svdof, vdof2 in zip(svdofs, vdofs2):
                x_merge[svdof] = gf[vdof2]
            iel = iel+1
    totals = x_merge.GetDataArray()
    print("merged gridfunction in myid=0", totals)

    x_merge.Save("merge_sol_merged.gf")

    m = mfem.Mesh('merge_sol.mesh', 1, 0, False)
    gf_ref = mfem.GridFunction(m, 'merge_sol_ref.gf')
    print(gf_ref.GetDataArray())
    print("differnece", np.sum(np.abs(gf_ref.GetDataArray() - totals)))

def repartition_gf(parts):
    mesh = mfem.Mesh('merge_sol.mesh', 1, 0, False)
    gf = mfem.GridFunction(mesh, 'merge_sol_merged.gf')

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh, parts)
    fec = fec_type(fec_order, mesh.Dimension())        
    pfes = mfem.ParFiniteElementSpace(pmesh, fec, 1)        

    pgf = mfem.ParGridFunction(pmesh, gf, parts)
    smyid = '.{:0>6d}'.format(myid)            
    pgf.Save('merge_sol_data2'+smyid, 16)

def run_test():
    # Generate sample data
    parts = generate_data()

    # Merge data on root node.
    if myid == 0:
        merge_data(parts)
    MPI.COMM_WORLD.Barrier()
    repartition_gf(parts)
    
if __name__=='__main__':
    run_test()
