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
    sigma_attr_coefs = mfem.MatrixCoefficientPtrArray(max_attr)
    sigma_attr = mfem.intArray(max_attr)
    tensors =  [mfem.DenseMatrix(np.ones((3,3))*i) for i in range(max_attr)]
    tensor_coefs =  [mfem.MatrixConstantCoefficient(mat) for mat in tensors]
    
    for ti, tensor in enumerate(tensors): 
        # add matrix coefficient to list
        print(ti)
        xx = mfem.MatrixConstantCoefficient(tensor)
        sigma_attr_coefs[ti] = xx
        sigma_attr[ti] = ti+1   
    
    # Create PW Matrix Coefficient
    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, sigma_attr_coefs, False)
    sigmaCoef = mfem.PWMatrixCoefficient(dim, sigma_attr, tensor_coefs, False)


if __name__ == '__main__':
    test()
