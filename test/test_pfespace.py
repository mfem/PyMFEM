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
    from mfem.common.parcsr_extra import ToScipyCSR, ToScipyCoo, ReadHypreDiag
else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0


def get_inv_doftrans(doftrans):
    if doftrans is not None:
        Mdoftrans = np.zeros((doftrans.Size(), doftrans.Size()))
        vv = mfem.Vector(doftrans.Size())
        for i in range(doftrans.Size()):
            vv.Assign(0)
            vv[i] = 1
            doftrans.InvTransformPrimal(vv)
            Mdoftrans[:, i] = vv.GetDataArray()

        return Mdoftrans
    else:
        return None


def run_test():
    order = 1
    vdim = 1
    path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    meshfile = os.path.join(path, 'data', 'waveguide.mesh')

    mesh = mfem.Mesh(meshfile, 1, 1)

    dim = mesh.Dimension()
    sdim = mesh.SpaceDimension()
    fec = mfem.ND_FECollection(order, sdim)
    #fec = mfem.H1_FECollection(order, sdim)

    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    # pmesh.ReorientTetMesh()
    pfes = mfem.ParFiniteElementSpace(pmesh, fec, vdim)
    print('Number of global finite element unknowns: ' +
          str(pfes.GetTrueVSize()))

    P = pfes.Dof_TrueDof_Matrix()

    print("Is P HypreGPU?s:", mfem.is_HYPRE_USING_CUDA())

    print("Height/Width", P.Height(), P.Width())

    cparts = P.GetColPartArray()
    rparts = P.GetRowPartArray()
    m, n = P.M(), P.N()

    print("col part", cparts)
    print("row part", rparts)

    # this does not produce correct diag for non-square matrix !?
    #diag = mfem.Vector()
    # P.GetDiag(diag)
    #d1 = diag.GetDataArray()

    #d2 = ReadHypreDiag(P)
    #print("diag length:", d2.shape[0])

    merged = mfem.SparseMatrix()
    P.MergeDiagAndOffd(merged)

    print("checking number of rows: ", merged.Height())
    assert merged.Height() == P.Height(), "local height does not match"
    print("passed")
    #print("I", merged.GetIArray())
    #print("J", merged.GetJArray())

    print("checking length of local nnz: ", P.get_local_nnz())
    assert len(merged.GetDataArray()
               ) == P.get_local_nnz(), "local nnz not correct"
    print("passed")

    print("checking NNZ: ", P.NNZ())
    assert np.sum(MPI.COMM_WORLD.allgather(P.get_local_nnz())
                  ) == P.NNZ(), "total NNZ does not match"
    print("passed")

    from mfem.common.sparse_utils import sparsemat_to_scipycsr

    merged_csr = sparsemat_to_scipycsr(merged, float)

    # if mfem.is_HYPRE_USING_CUDA():
    Pcsr = ToScipyCSR(P)
    Pcsr2 = ToScipyCoo(P).tocsr()

    print("checking various matrix conversion")
    assert (merged_csr != Pcsr).nnz == 0, "ToScipyCSR deos not agree with merged matrix"
    assert (merged_csr !=
            Pcsr2).nnz == 0, "ToScipyCoo deos not agree with merged matrix"
    print("passed")

    print("size", Pcsr.indptr.shape, Pcsr.shape)
    VDoFtoGTDoF1 = Pcsr.indices

    print("Global V and True sizes", pfes.GlobalVSize(), pfes.GlobalTrueVSize())

    # print(pfes.GetNE())
    # for i in range(pfes.GetNE()):
    for i in [1]:
        vdofs = pfes.GetElementDofs(i)
        doftrans = pfes.GetElementDofTransformation(i)
        vdofs = [-1 - x if x < 0 else x for x in vdofs]
        ltdof = [pfes.GetLocalTDofNumber(ii) for ii in vdofs]
        gtdof1 = [list(Pcsr.indices[Pcsr.indptr[ii:ii+1]])
                  for ii in vdofs]  # if pfes.GetLocalTDofNumber(ii) >= 0]
        gtdof1 = sum(gtdof1, [])
        # if pfes.GetLocalTDofNumber(ii) >= 0]
        gtdof2 = [pfes.GetGlobalTDofNumber(ii) for ii in vdofs]
        print(gtdof1, gtdof2)
        print(gtdof1 != gtdof2)
        if True:
            print("myid: ", myid)
            # print(get_inv_doftrans(doftrans))
            print("Local DoFs", vdofs)
            print("Local TDoFs", ltdof)
            print("Global TDoFs from P             :", gtdof1)
            print("Global TDoFs GetGlobalTDoFNumber:", gtdof2)


if __name__ == "__main__":

    if not use_parallel:
        print("this is parallel only, skipping")
    else:
        run_test()
