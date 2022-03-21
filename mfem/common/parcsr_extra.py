from __future__ import print_function

#
#   this modules works only with parallel version
#
import sys
import numpy as np


def get_assumed_patitioning(m):
    '''
    for given size of row, returns proper patitioning
    '''
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    min_nrows = m // num_proc
    extra_rows = m % num_proc
    start_row = min_nrows * myid + (extra_rows if extra_rows < myid else myid)
    end_row = start_row + min_nrows + (1 if extra_rows > myid else 0)
    nrows = end_row - start_row

    return start_row, end_row, nrows


def get_row_partitioning(M):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    m = M.GetNumRows()
    m_array = comm.allgather(m)
    rows = [0] + list(np.cumsum(m_array))
    return rows


def get_col_partitioning(M):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    m = M.GetNumCols()
    m_array = comm.allgather(m)
    rows = [0] + list(np.cumsum(m_array))
    return rows


def ToHypreParVec(vec):
    import mfem.par as mfem
    from mpi4py import MPI

    if mfem.sizeof_HYPRE_Int() == 4:
        dtype = 'int32'
    else:
        dtype = 'int64'

    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    vec = vec.flatten()
    ml = vec.shape[0]

    # collect col array to determin partitioning
    m_array = comm.allgather(ml)
    cols = [0] + list(np.cumsum(m_array))
    glob_size = cols[-1]
    col_starts = np.array([cols[myid], cols[myid+1], glob_size], dtype=dtype)

    vec = vec.astype('float', copy=False)
    v = mfem.HypreParVector(MPI.COMM_WORLD,
                            glob_size, [vec, col_starts])

    return v


def ToHypreParCSR(mat, check_partitioning=False, verbose=False,
                  col_starts=None, assert_non_square_no_col_starts=True):
    '''
    convert scipy sparse matrix to hypre

    vertically stack csr matrix to generte HYPRE Par CSR

    Note:
    row partitioning is inferred from distribution of input matrix.
    column patitioning needs to be specified col_starts.

    If col_starts is not given, column partitioning is chosen 
    to be the same as row partitioning. This works if matrix is square (M = N).

    For an aribtrary rectangular matrix, the column partitioning can be
    different from the row partitioning. For example, MFEM mixedbilinearfomr 
    uses different partitiong rules for row and column.

    ToDo: change default assert_non_square_no_col_starts to False

    '''

    from mpi4py import MPI
    import mfem.par as mfem

    if mfem.sizeof_HYPRE_Int() == 4:
        dtype = 'int32'
    else:
        dtype = 'int64'

    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    def verbose_message(m, n, nrows, i, j, data, row_starts, col_starts):
        for k in range(num_proc):
            MPI.COMM_WORLD.Barrier()
            if myid == k:
                print('MyID : ', k)
                print((m, n), nrows, len(data), i, j,
                      data, row_starts, col_starts)
                print('NNZ', np.sum(data != 0.0))
        MPI.COMM_WORLD.Barrier()

    from scipy.sparse import csr_matrix

    if isinstance(mat, csr_matrix):
        mat = mat.astype('float')
        ml, nl = mat.shape
        n_array = comm.allgather(nl)
    else:
        raise ValueError("Import Matrix Format should be csr or None")

    # collect row array to determin the size of matrix
    m_array = comm.allgather(ml)

    rows = [0] + list(np.cumsum(m_array))
    m = rows[-1]
    #row_starts = np.array([rows[myid], rows[myid+1], m], dtype=dtype)
    row_starts = np.array([rows[myid], rows[myid+1]], dtype=dtype)

    n = nl
    nrows = ml

    i = mat.indptr.astype(dtype)
    j = mat.indices.astype(dtype)
    data = mat.data

    if col_starts is None and m != nl:
        col_starts = get_assumed_patitioning(nl)
        if assert_non_square_no_col_starts:
            assert False, "col_starts must be specified for non diagonal array"
    if col_starts is None:
        col_starts = row_starts.copy()
        # col_starts[-1]=n
        if col_starts[0] > n:
            col_starts[0] = n
        if col_starts[1] > n:
            col_starts[1] = n
        #col_starts[2] = n
    else:
        # make sure that dtype is right....
        col_starts = np.array(col_starts, dtype=dtype)
    if check_partitioning:
        ch = get_assumed_patitioning(m)
        if (row_starts[0] != ch[0] or
            row_starts[1] != ch[1] or
                nrows != ch[2]):
            for k in range(num_proc):
                MPI.COMM_WORLD.Barrier()
                if myid == k:
                    print('MyID : ', k)
                    print(ch, nrows, row_starts, col_starts)
                    print('NNZ', np.sum(data != 0.0))
            MPI.COMM_WORLD.Barrier()
            raise ValueError("partitioning of input matrix is not correct")
    if verbose:
        verbose_message(m, n, nrows, i, j, data, row_starts, col_starts)

    #
    # it seems row_starts and col_starts are both to determin
    # which part is treated diagnal element.
    #
    if (m == n and row_starts[0] == col_starts[0] and
            row_starts[1] == col_starts[1]):
        # this will cause hypre_CSRMatrixReorder call.
        M = mfem.HypreParMatrix(MPI.COMM_WORLD,
                                nrows,
                                m, n, [i, j,
                                       data, col_starts])
        M.CopyRowStarts()
        M.CopyColStarts()
    else:
        M = mfem.HypreParMatrix(MPI.COMM_WORLD,
                                nrows,
                                m, n, [i, j,
                                       data, row_starts[:2], col_starts[:2]])
        M.CopyRowStarts()
        M.CopyColStarts()
    return M


def ToScipyCSR(mat):
    '''
    convert HypreParCSR to Scipy CSR Matrix
    '''
    import mfem.par as mfem
    merged = mfem.SparseMatrix()
    mat.MergeDiagAndOffd(merged)

    from scipy.sparse import csr_matrix
    P = csr_matrix((merged.GetDataArray(), merged.GetJArray(), merged.GetIArray()),
                   shape=(merged.Height(), merged.Width()))
    P._linked_mat = merged
    return P


def ToScipyCoo(mat):
    '''
    convert HypreParCSR to Scipy Coo Matrix
    '''
    num_rows, ilower, iupper, jlower, jupper, irn, jcn, data = mat.GetCooDataArray()
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    m = iupper - ilower + 1
    n = jupper - jlower + 1
    n = mat.N()

    from scipy.sparse import coo_matrix
    try:
        return coo_matrix((data, (irn-ilower, jcn)), shape=(m, n))
    except:
        print("wrong input")
        print(num_rows, ilower, iupper, jlower, jupper)
        print(np.min(irn-ilower), np.max(irn-ilower),
              np.min(jcn),  np.max(jcn), (m, n))
        raise


def InnerProductComplex(A, B):
    def ensure_hypreread(V):
        if V is None:
            return
        if not V._hypreread_called:
            V.HypreRead()

    import mfem.par as mfem
    R_A, I_A = A
    R_B, I_B = B

    for V in (R_A, I_A, R_B, I_B):
        ensure_hypreread(V)

    if I_A is None and I_B is None:
        return mfem.InnerProduct(R_A, R_B)
    elif I_A is None:
        r = mfem.InnerProduct(R_A, R_B)
        i = mfem.InnerProduct(R_A, I_B)
    elif I_B is None:
        r = mfem.InnerProduct(R_A, R_B)
        i = mfem.InnerProduct(I_A, R_B)
    else:
        r = mfem.InnerProduct(R_A, R_B) - mfem.InnerProduct(I_A, I_B)
        i = mfem.InnerProduct(R_A, I_B) + mfem.InnerProduct(I_A, R_B)
    return r + 1j * i


def ParAdd(A, B):
    '''
    add HypreParCSR

    '''
    col_starts = A.GetColPartArray()  # ; col_starts[2] = A.N()
    return ToHypreParCSR((ToScipyCoo(A) + ToScipyCoo(B)).tocsr(),
                         col_starts=col_starts)


def ParMultVecComplex(A, v):
    '''
    A*v
    '''
    import mfem.par as mfem
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    R_A, I_A = A
    R_v, I_v = v

    # Take Row partitioning of A for output
    if R_A is not None:
        part = R_A.GetRowPartArray()
    elif I_A is not None:
        part = I_A.GetRowPartArray()
    else:
        return (None, None)

    ans_r = ToHypreParVec(np.zeros(part[1]-part[0]))
    if I_A is None and I_v is None:
        R_A.Mult(R_v, ans_r)
        return (ans_r, None)
    else:
        ans_i = ToHypreParVec(np.zeros(part[1]-part[0]))

    if I_A is None:
        R_A.Mult(R_v, ans_r)
        R_A.Mult(I_v, ans_i)
    elif I_v is None:
        R_A.Mult(R_v, ans_r)
        I_A.Mult(R_v, ans_i)
    else:
        ans_r2 = ToHypreParVec(np.zeros(part[1]-part[0]))
        ans_i2 = ToHypreParVec(np.zeros(part[1]-part[0]))
        R_A.Mult(R_v, ans_r)
        I_A.Mult(I_v, ans_r2)
        ans_r -= ans_r2

        R_A.Mult(I_v, ans_i)
        I_A.Mult(R_v, ans_i2)
        ans_i += ans_i2

    return (ans_r, ans_i)


def ParMultComplex(A, B):
    '''
    compute complex mult of hypre real matrices

    A = (R_A, I_A)
    B = (R_B, I_B)

    (R_A*R_B - I_A*I_B, R_A*I_B + I_A*R_B)
    '''
    from mpi4py import MPI
    import mfem.par as mfem

    comm = MPI.COMM_WORLD
    num_proc = MPI.COMM_WORLD.size
    myid = MPI.COMM_WORLD.rank

    R_A, I_A = A
    R_B, I_B = B

    if I_A is None and I_B is None:
        r = mfem.ParMult(R_A, R_B)
        r.CopyRowStarts()
        r.CopyColStarts()

        return (r, None)
    elif I_A is None:
        r = mfem.ParMult(R_A, R_B)
        i = mfem.ParMult(R_A, I_B)
        r.CopyRowStarts()
        r.CopyColStarts()
        i.CopyRowStarts()
        i.CopyColStarts()
        return (r, i)

    elif I_B is None:
        r = mfem.ParMult(R_A, R_B)
        i = mfem.ParMult(I_A, R_B)
        r.CopyRowStarts()
        r.CopyColStarts()
        i.CopyRowStarts()
        i.CopyColStarts()

        return (r, i)
    else:
        A = mfem.ParMult(R_A, R_B)
        B = mfem.ParMult(I_A, I_B)
        C = mfem.ParMult(R_A, I_B)
        D = mfem.ParMult(I_A, R_B)
        col_starts = A.GetColPartArray()  # ; col_starts[2] = A.N()
        r = ToHypreParCSR((ToScipyCoo(A) - ToScipyCoo(B)).tocsr(),
                          col_starts=col_starts)
        i = ToHypreParCSR((ToScipyCoo(C) + ToScipyCoo(D)).tocsr(),
                          col_starts=col_starts)
        return (r, i)


def TransposeComplex(A):
    '''
    A is tuple (A_real, A_imag), whose real/imag are
    HypreParCSR
    '''

    R = A[0].Transpose() if A[0] is not None else None
    I = A[1].Transpose() if A[1] is not None else None
    if R is not None:
        R.thisown = True
    if I is not None:
        I.thisown = True
    return (R, I)


def Conj(A):

    R = A[0]
    I = A[1]
    if I is None:
        return A
    col_starts = I.GetColPartArray()
    col_starts[2] = I.N()
    I = ToHypreParCSR(-ToScipyCoo(I).tocsr(), col_starts=col_starts)
    return (R, I)


def RapComplex(A, B):
    '''
    Bt * A * B

    for complex A and B
    '''
    X = ParMultComplex(A, B)
    return ParMultComplex(Conj(TransposeComplex(B)), X)


def Array2HypreVec(v, partitioning=None, rank=0):
    '''
    convert array in rank (default = 0)  to 
    distributed Hypre 1D Matrix (size = m x 1)
    '''
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()

    data = v if myid == rank else None
    data = comm.bcast(data, root=rank)

    if partitioning is None:
        start_row, end_row, nrows = get_assumed_patitioning(len(data))
    else:
        start_row = partitioning[myid]
        end_row = partitioning[myid+1]
        nrows = end_row - start_row

    from scipy.sparse import csr_matrix, coo_matrix
    v = np.ascontiguousarray(data[start_row:end_row].flatten())
    return ToHypreParVec(v)
    #m = csr_matrix(np.array(d).reshape(-1,1), shape=(nrows,1), dtype='float')
    # return ToHypreParCSR(m)


def HypreVec2Array(V, copy=True):
    '''
    convert HypreParVec to 1D array 
    on rank = 0
    '''
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    data = V.GetDataArray()
    if copy:
        data = data.copy()

    rcounts = len(data)

    rcounts = MPI.COMM_WORLD.gather(rcounts, root=0)
    cm = np.hstack((0, np.cumsum(rcounts)))
    disps = list(cm[:-1])
    recvdata = None
    senddata = [data, data.shape[0]]

    if myid == 0:
        length = cm[-1]
        recvbuf = np.empty([length], dtype='float')
        recvdata = [recvbuf, rcounts, disps, MPI.DOUBLE]
    else:
        recvdata = [None, rcounts, disps, MPI.DOUBLE]
        recvbuf = None

    MPI.COMM_WORLD.Gatherv(senddata, recvdata,  root=0)
    if myid == 0:
        MPI.COMM_WORLD.Barrier()
        return np.array(recvbuf)
    MPI.COMM_WORLD.Barrier()
    return None


def ResetHypreDiag(M, idx, value=1.0):
    '''
    set diagonal element to value (normally 1)
    '''
    col_starts = M.GetColPartArray()

    num_rows, ilower, iupper, jlower, jupper, irn, jcn, data = M.GetCooDataArray()
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    m = iupper - ilower + 1
    n = jupper - jlower + 1
    n = M.N()
    from scipy.sparse import coo_matrix, lil_matrix
    try:
        mat = coo_matrix((data, (irn-ilower, jcn)), shape=(m, n)).tolil()
    except:
        print("wrong input")
        print(num_rows, ilower, iupper, jlower, jupper)
        print(np.min(irn-ilower), np.max(irn-ilower),
              np.min(jcn),  np.max(jcn), (m, n))
        raise

    idx = np.array(idx, dtype=int, copy=False)
    ii = idx[np.logical_and(idx >= ilower, idx <= iupper)]
    mat[ii-ilower, ii] = value
    # for ii in idx:
    #    if ii >= ilower and ii <= iupper:
    #       mat[ii-ilower, ii] = value

    return ToHypreParCSR(mat.tocsr(), col_starts=col_starts)


def ResetHypreRow(M, idx):
    '''
    set row 0.0
    '''
    col_starts = M.GetColPartArray()  # ; col_starts[2] = M.N()
    num_rows, ilower, iupper, jlower, jupper, irn, jcn, data = M.GetCooDataArray()

    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    m = iupper - ilower + 1
    n = jupper - jlower + 1
    n = M.N()
    from scipy.sparse import coo_matrix, lil_matrix

    k = np.in1d(irn, idx)
    data[k] = 0.0

    mat = coo_matrix((data, (irn-ilower, jcn)), shape=(m, n)).tocsr()
    mat.eliminate_zeros()

    return ToHypreParCSR(mat.tocsr(), col_starts=col_starts)


def ResetHypreCol(M, idx):
    '''
    set col zero
    '''
    col_starts = M.GetColPartArray()  # ; col_starts[2] = M.N()
    num_rows, ilower, iupper, jlower, jupper, irn, jcn, data = M.GetCooDataArray()

    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    m = iupper - ilower + 1
    n0 = jupper - jlower + 1
    n = M.N()
    from scipy.sparse import coo_matrix, lil_matrix

    k = np.in1d(jcn, idx)
    data[k] = 0.0

    mat = coo_matrix((data, (irn-ilower, jcn)), shape=(m, n)).tocsr()
    mat.eliminate_zeros()
    return ToHypreParCSR(mat.tocsr(), col_starts=col_starts)


def ReadHypreDiag(M, idx):
    '''
    get diagonal element
    '''
    col_starts = M.GetColPartArray()

    num_rows, ilower, iupper, jlower, jupper, irn, jcn, data = M.GetCooDataArray()
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    m = iupper - ilower + 1
    n = M.N()
    from scipy.sparse import coo_matrix, lil_matrix
    try:
        mat = coo_matrix((data, (irn-ilower, jcn)), shape=(m, n)).tolil()
    except:
        print("wrong input")
        print(num_rows, ilower, iupper, jlower, jupper)
        print(np.min(irn-ilower), np.max(irn-ilower),
              np.min(jcn),  np.max(jcn), (m, n))
        raise

    idx = np.array(idx, dtype=int, copy=False)
    ii = idx[np.logical_and(idx >= ilower, idx <= iupper)]

    tmp = mat[ii-ilower, ii].toarray().flatten()

    #idx = np.arange(ilower, min([iupper+1, n]))
    #ii = idx[np.logical_and(idx >= ilower, idx <= iupper)]
    #tmp = mat[ii-ilower, ii].toarray().flatten()

    return tmp
