import sys
import numpy as np
import gc
from scipy.sparse import csr_matrix

if len(sys.argv) > 1 and sys.argv[1] == '-p':
    import mfem.par as mfem
    use_parallel = True

    from mfem.common.mpi_debug import nicePrint as print

    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0

def test():
    dim = 3
    max_attr = 5
    sigma_attr_coefs = mfem.MatrixCoefficientArray(max_attr)
    sigma_attr = mfem.intArray(max_attr)


    tensors =  [mfem.DenseMatrix(np.ones((3,3))*i) for i in range(max_attr)]
    tensor_coefs =  [mfem.MatrixConstantCoefficient(mat) for mat in tensors]
    sigma_array = mfem.MatrixCoefficientArray(tensor_coefs)
    sigma_attr = mfem.intArray([1,2,3,4,5])
    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, sigma_array, False)

    for ti, tensor in enumerate(tensors): 
        # add matrix coefficient to list
        xx = mfem.MatrixConstantCoefficient(tensor)
        sigma_attr_coefs[ti] = xx
        sigma_attr[ti] = ti+1   

    # Create PW Matrix Coefficient
    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, sigma_attr_coefs, False)

    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, tensor_coefs, False)
    tensor_coefs = mfem.MatrixCoefficientArray([mfem.MatrixConstantCoefficient(mat) for mat in tensors])
    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, tensor_coefs, False)

    data = tensor_coefs.GetData()
    tensor_coefs2 = mfem.MatrixCoefficientArray(data, 5, False)
    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, tensor_coefs2, False)

    print("exiting")


if __name__ == '__main__':
    import tracemalloc

    tracemalloc.start()

    print(tracemalloc.get_traced_memory())

    for i in range(3000):
        test()
        print(tracemalloc.get_traced_memory())
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')
    for stat in top_stats[:10]: 
        print(stat)

