from __future__ import print_function
import os
import sys
import numpy as np

from scipy.sparse import csr_matrix, coo_matrix


def run_test():
    import mfem.par as par
    from mfem.common.parcsr_extra import ToHypreParCSR, ToScipyCoo
    from mpi4py import MPI
    from mfem.common.mpi_debug import nicePrint

    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    device = mfem.Device('cuda')
    
    def print_hypre(M, txt):
        for i in range(num_proc):
            MPI.COMM_WORLD.Barrier()
            if myid == i:
                if myid == 0:
                    print(txt)
                    print('MyID: ', myid)
                else:
                    print('MyID: ', myid)
                print(ToScipyCoo(M))

    # make sample matrix
    row = np.array([0, 0, 1, 1])
    col = np.array([0, 3, 1, 2])
    data = np.array([4, 5, 7, 9])
    m = coo_matrix((data, (row, col)), shape=(2, 4))
    m = m.tocsr()
    m = m*(myid+1)

    M = ToHypreParCSR(m, assert_non_square_no_col_starts=False, verbose=True)

    print("memory class, memory type",
          mfem.GetHypreMemoryClass(), mfem.GetHypreMemoryType())

    print_hypre(M, 'matrix M')

    from mfem.common.chypre import CHypreVec

    r1 = np.array([0, 0, 1, 1])
    r2 = np.array([1, 1, 0, 0])
    vec1 = CHypreVec(r1, None)
    vec2 = CHypreVec(r2, None)

    if myid == 0:
        print("v1")
    v1 = (vec1-vec1*1j)
    v2 = (vec1+vec1*1j)

    nicePrint(v1.GlobalVector())
    nicePrint(v2.GlobalVector())
    nicePrint((v1+v2).GlobalVector())
    nicePrint((v1-v2).GlobalVector())

    if myid == 0:
        print("v1, v2")
    v1 = (vec1-vec2*1j)
    v2 = (vec1+vec2*1j)
    nicePrint(v1.GlobalVector())
    nicePrint(v2.GlobalVector())

    if myid == 0:
        print("v1+v2")
    nicePrint((v1+v2).GlobalVector())
    if myid == 0:
        print("v1-v2")
    nicePrint((v1-v2).GlobalVector())

    if myid == 0:
        print("3*v1")
    nicePrint(v1.GlobalVector())
    v1 *= 1j
    if myid == 0:
        print("1j*v1")
    nicePrint(v1.GlobalVector())
    if myid == 0:
        print("v1 dot v1")
    print(v1.dot(v1))
    if myid == 0:
        print("v1 *= 1+1j")
    v1 *= 1+1j
    nicePrint("v1", v1.GlobalVector())
    nicePrint("v2", v2.GlobalVector())
    #print(v1.dot(v1))
    #print(v1.dot(v2))
    #if mfem.is_HYPRE_USING_CUDA():
    #    v1.HypreRead()

    v1[0] += 3.
    if myid == 0:
        print("3*v1")
    #if mfem.is_HYPRE_USING_CUDA():
    #    v1.HypreRead()
        
    nicePrint(v1[0].GetDataArray())
    nicePrint(v1.GlobalVector())

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':
        import mfem.par as mfem
        run_test()
    else:
        import mfem.ser as mfem
