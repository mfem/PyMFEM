'''
  CHypre (Complex Hypre)

  CHypreVec : ParVector
  CHypreMat : ParCSR

  container object to support complex using
  real value hypre

  it should work with pure real or pure imaginary
  case too.

  it follows the mathod naming convetion used
  in scipy.sparse. However, since it inherits the list
  object, __setitem__ can not be used for accessing
  array  elements. Use set_element, instead.
'''
import numpy as np
from numbers import Number
from scipy.sparse import csr_matrix, coo_matrix, lil_matrix
from mfem.common.parcsr_extra import *

# DO NOT IMPORT MPI in Global, sicne some routins will be used
# in serial mode too.
try:
    import mfem.par
    MFEM_PAR = True
except BaseException:
    MFEM_PAR = False


class CHypreVec(list):
    def __init__(self, r=None, i=None, horizontal=False):
        list.__init__(self, [None] * 2)
        self._horizontal = horizontal

        if isinstance(r, np.ndarray):
            self[0] = ToHypreParVec(r)
        else:
            self[0] = r
        if isinstance(i, np.ndarray):
            self[1] = ToHypreParVec(i)
        else:
            self[1] = i

    def __repr__(self):
        if self[0] is not None:
            part = self[0].GetPartitioningArray()
        elif self[1] is not None:
            part = self[1].GetPartitioningArray()
        else:
            return "CHypreVec (empty)"
        return "CHypreVec" + str(self.shape) + \
            "[" + str(part[1] - part[0]) + "]"

    @property
    def imag(self):
        return self[1]

    @imag.setter
    def imag(self, value):
        self[1] = value

    @property
    def real(self):
        return self[0]

    @real.setter
    def real(self, value):
        self[0] = value

    @property
    def shape(self):
        if self[0] is not None:
            size = self[0].GlobalSize()
        elif self[1] is not None:
            size = self[1].GlobalSize()
        else:
            size = 0.0

        if self._horizontal:
            return 1, size
        else:
            return size, 1

    def HypreRead(self):
        if self[0] is not None:
            self[0].HypreRead()
        if self[1] is not None:
            self[1].HypreRead()

    def HypreReadWrite(self):
        if self[0] is not None:
            self[0].HypreReadWrite()
        if self[1] is not None:
            self[1].HypreReadWrite()

    def HypreWrite(self):
        if self[0] is not None:
            self[0].HypreWrite()
        if self[1] is not None:
            self[1].HypreWrite()

    def isComplex(self):
        return not (self[1] is None)

    def GetPartitioningArray(self):
        if self[0] is not None:
            part = self[0].GetPartitioningArray()
            #part[2] = self[0].GlobalSize()
        elif self[1] is not None:
            prat = self[1].GetPartitioningArray()
            #part[2] = self[1].GlobalSize()
        else:
            raise ValueError("CHypreVec is empty")
        return part

    def __imul__(self, other):
        if isinstance(other, CHypreVec):
            assert False, "CHypreVec *= vector is not supported. Use dot"
        elif np.iscomplexobj(other):
            #other = complex(other)
            i = other.imag
            r = other.real
            if self[0] is not None and self[1] is not None:
                rr = self[0].GetDataArray() * r - self[1].GetDataArray() * i
                ii = self[0].GetDataArray() * i + self[1].GetDataArray() * r
                self[0] = ToHypreParVec(rr)
                self[1] = ToHypreParVec(ii)

            elif self[0] is not None:
                if np.any(i != 0.):
                    self[1] = ToHypreParVec(i * self[0].GetDataArray())
                if np.any(r != 0.):
                    tmp = self[0].GetDataArray()
                    tmp *= r
                else:
                    self[0] = None

            elif self[1] is not None:
                if np.any(i != 0.):
                    self[0] = ToHypreParVec(-i * self[1].GetDataArray())
                if np.any(r != 0.):
                    tmp = self[1].GetDataArray()
                    tmp *= r
                else:
                    self[1] = None

            else:
                passself[0] = None
        else:
            other = float(other)
            if self[0] is not None:
                self[0] *= other
            if self[1] is not None:
                self[1] *= other
        return self

    def __mul__(self, other):
        if isinstance(other, CHypreVec):
            assert False, "CHypreVec *= vector is not supported. Use dot"
        elif np.iscomplexobj(other):
            other = complex(other)
            i = other.imag
            r = other.real
        else:
            r = float(other)
            i = 0.0
        rdata = self[0].GetDataArray() if self[0] is not None else 0
        idata = self[1].GetDataArray() if self[1] is not None else 0

        rr = rdata * r - idata * i
        ii = rdata * i + idata * r

        # note: for the real part we keep it even if it is zero
        #       so that it conservs vector size information
        rr = ToHypreParVec(rr)
        ii = ToHypreParVec(ii) if np.count_nonzero(ii) != 0 else None

        return CHypreVec(rr, ii, horizontal=self._horizontal)

    def __add__(self, other):
        assert self._horizontal == other._horizontal, "can not add vertical and hirozontal vector"

        if self[0] is not None and other[0] is not None:
            data = self[0].GetDataArray() + other[0].GetDataArray()
            r = ToHypreParVec(data)
        elif self[0] is not None:
            data = self[0].GetDataArray()
            r = ToHypreParVec(data)
        elif other[0] is not None:
            data = other[0].GetDataArray()
            r = ToHypreParVec(data)
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            data = self[1].GetDataArray() + other[1].GetDataArray()
            i = ToHypreParVec(data)
        elif self[1] is not None:
            data = self[1].GetDataArray()
            i = ToHypreParVec(data)
        elif other[1] is not None:
            data = other[1].GetDataArray()
            i = ToHypreParVec(data)
        else:
            i = None
        return CHypreVec(r, i, horizontal=self._horizontal)

    def __sub__(self, other):
        assert self._horizontal == other._horizontal, "can not add vertical and hirozontal vector"
        if self[0] is not None and other[0] is not None:
            data = self[0].GetDataArray() - other[0].GetDataArray()
            r = ToHypreParVec(data)
        elif self[0] is not None:
            data = self[0].GetDataArray()
            r = ToHypreParVec(data)
        elif other[0] is not None:
            data = -other[0].GetDataArray()
            r = ToHypreParVec(data)
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            data = self[1].GetDataArray() - other[1].GetDataArray()
            i = ToHypreParVec(data)
        elif self[1] is not None:
            data = self[1].GetDataArray()
            i = ToHypreParVec(data)
        elif other[1] is not None:
            data = -other[1].GetDataArray()
            i = ToHypreParVec(data)
        else:
            i = None
        return CHypreVec(r, i, horizontal=self._horizontal)

    def dot(self, other):
        if isinstance(other, CHypreVec):
            return InnerProductComplex(self, other)
        elif (isinstance(other, CHypreMat) and
              self._horizontal):
            ret = other.transpose().dot(self)
            ret._horizontal = True
            return ret
        else:
            raise ValueError(
                "CHypreVec::dot supports Vec*Vec (InnerProduct) and (Mat^t*Vec)^t ")

    def get_elements(self, idx):
        part = self.GetPartitioningArray()
        idx = idx - part[0]
        idx = idx[idx < part[1]-part[0]]
        idx = idx[idx >= 0]

        if len(idx) == 0:
            return np.array([])

        ret = 0.0
        if self[0] is not None:
            ret = ret + self[0].GetDataArray()[idx]
        if self[1] is not None:
            ret = ret + 1j*self[1].GetDataArray()[idx]
        return ret

    def set_elements(self, idx, value):
        part = self.GetPartitioningArray()
        idx = idx - part[0]
        idx = idx[idx < part[1]-part[0]]
        idx = idx[idx >= 0]

        rvalue = value.real if np.iscomplexobj(value) else value

        if len(idx) == 0:
            return

        if self[0] is not None:
            self[0].GetDataArray()[idx] = rvalue

        if np.iscomplexobj(value):
            if self[1] is None:
                i = self[0].GetDataArray()*0.0
                self[1] = ToHypreParVec(i)
            self[1].GetDataArray()[idx] = value.imag

    def set_element(self, i, v):
        part = self.GetPartitioningArray()
        if part[0] <= i and i < part[1]:
            v = complex(v)
            if self[0] is not None:
                self[0][int(i - part[0])] = v.real
            if self[1] is not None:
                self[1][int(i - part[0])] = v.imag

    def get_element(self, i):
        part = self.GetPartitioningArray()
        if part[0] <= i and i < part[1]:
            if self[0] is not None:
                r = self[0][int(i - part[0])]
            else:
                r = 0
            if self[1] is not None:
                return r + 1j * self[1][int(i - part[0])]
            else:
                return r

    def copy_element(self, tdof, vec):
        for i in tdof:
            v = vec.get_element(i)
            self.set_element(i, v)

    '''
    def gather(self):
        from mpi4py import MPI
        myid = MPI.COMM_WORLD.rank
        vecr = 0.0; veci = 0.0
        if self[0] is not None:
           vecr  = gather_vector(self[0].GetDataArray(), MPI.DOUBLE)
        if self[1] is not None:
           veci = gather_vector(self[1].GetDataArray(), MPI.DOUBLE)

        if myid == 0:
            if self[0] is None:
                return vecr
            else:
                return vecr + 1j*veci
    '''

    def get_squaremat_from_right(self):
        '''
        squre matrix which can be multipled from the right of self.
        '''
        if not self._horizontal:
            raise ValueError("Vector orientation is not right")

        part = self.GetPartitioningArray()
        width = self[1].GlobalSize()
        return SquareCHypreMat(width, part, real=(self[1] is None))

    def transpose(self):
        self._horizontal = not self._horizontal
        return self

    def _do_reset(self, v, idx):
        # ownership is transferrd to a new vector.
        ownership = v.GetOwnership()
        data = v.GetDataArray()
        part = v.GetPartitioningArray()
        for i in idx:
            if i >= part[1]:
                continue
            if i < part[0]:
                continue
            data[i - part[0]] = 0

        ret = ToHypreParVec(data)
        ret.SetOwnership(ownership)
        v.SetOwnership(0)
        return ret

    def resetCol(self, idx):
        if self._horizontal:
            if self[0] is not None:
                self[0] = self._do_reset(self[0], idx)
            if self[1] is not None:
                self[1] = self._do_reset(self[1], idx)
        else:
            if 0 in idx:
                self *= 0.0

    def resetRow(self, idx):
        if self._horizontal:
            if 0 in idx:
                self *= 0.0
        else:
            if self[0] is not None:
                self[0] = self._do_reset(self[0], idx)
            if self[1] is not None:
                self[1] = self._do_reset(self[1], idx)

    def _do_select(self, v, lidx):
        # ownership is transferrd to a new vector.
        ownership = v.GetOwnership()
        data = v.GetDataArray()
        data2 = data[lidx]
        ret = ToHypreParVec(data2)
        ret.SetOwnership(ownership)
        v.SetOwnership(0)
        return ret

    def selectRows(self, idx):
        '''
        idx is global index
        '''
        if self._horizontal:
            if not 0 in idx:
                raise ValueError("VectorSize becomes zero")
            return self

        part = self.GetPartitioningArray()
        idx = idx[idx >= part[0]]
        idx = idx[idx < part[1]]
        lidx = idx - part[0]

        r = None
        i = None
        if not self._horizontal:
            if self[0] is not None:
                r = self._do_select(self[0], lidx)
            if self[1] is not None:
                i = self._do_select(self[1], lidx)
        return CHypreVec(r, i, horizontal=self._horizontal)

    def selectCols(self, idx):
        '''
        idx is global index
        '''
        if not self._horizontal:
            if not 0 in idx:
                raise ValueError("VectorSize becomes zero")
            return self

        part = self.GetPartitioningArray()
        idx = idx[idx >= part[0]]
        idx = idx[idx < part[1]]
        lidx = idx - part[0]

        r = None
        i = None
        if self._horizontal:
            if self[0] is not None:
                r = self._do_select(self[0], lidx)
            if self[1] is not None:
                i = self._do_select(self[1], lidx)
        return CHypreVec(r, i, horizontal=self._horizontal)

    def GlobalVector(self):
        '''
        Here it is reimplmente using MPI allgather.
        This is because GlobalVactor does not work when the vector
        is too small so that some node does not have data.
        '''
        def gather_allvector(data):
            from mfem.common.mpi_dtype import get_mpi_datatype
            from mpi4py import MPI

            mpi_data_type = get_mpi_datatype(data)
            myid = MPI.COMM_WORLD.rank
            rcounts = data.shape[0]
            rcounts = np.array(MPI.COMM_WORLD.allgather(rcounts))

            for x in data.shape[1:]:
                rcounts = rcounts * x
            cm = np.hstack((0, np.cumsum(rcounts)))
            disps = list(cm[:-1])
            length = cm[-1]
            recvbuf = np.empty([length], dtype=data.dtype)
            recvdata = [recvbuf, rcounts, disps, mpi_data_type]
            senddata = [data.flatten(), data.flatten().shape[0]]
            MPI.COMM_WORLD.Allgatherv(senddata, recvdata)
            return recvbuf.reshape(-1, *data.shape[1:])

        if self[0] is not None:
            v1 = gather_allvector(self[0].GetDataArray().copy())
        else:
            v1 = 0.0
        if self[1] is not None:
            v2 = gather_allvector(self[1].GetDataArray().copy())
            v1 = v1 + 1j * v2
        return v1

    def toarray(self):
        '''
        numpy array of local vector
        '''
        if self[0] is not None:
            v = self[0].GetDataArray()
        else:
            v = 0.0
        if self[1] is not None:
            v = v + 1j * self[1].GetDataArray()
        return v

    def isAllZero(self):
        return any(self.GlobalVector())

    def get_global_coo(self):
        data = self[0].GetDataArray()
        if self[1] is not None:
            data = data + 1j * self[1].GetDataArray()
        gcoo = coo_matrix(self.shape, dtype=data.dtype)
        gcoo.data = data
        part = self.GetPartitioningArray()
        if self._horizontal:
            gcoo.row = np.zeros(len(data))
            gcoo.col = np.arange(len(data)) + part[0]
        else:
            gcoo.col = np.zeros(len(data))
            gcoo.row = np.arange(len(data)) + part[0]
        return gcoo

    def true_nnz(self):
        '''
        more expensive version which reports nnz after
        eliminating all zero entries
        I can not use @property, since ScipyCoo can't
        have @property (class inherited from ndarray)
        '''
        data = self[0].GetDataArray()
        if self[1] is not None:
            data = data + 1j * self[1].GetDataArray()
        local_nnz = np.sum(data == 0.)

        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        nnz = comm.allgather(local_nnz)
        return nnz

    @property
    def isHypre(self):
        return True


class CHypreMat(list):
    def __init__(self, r=None, i=None, col_starts=None):
        list.__init__(self, [None] * 2)
        if isinstance(r, csr_matrix):
            self[0] = ToHypreParCSR(r, col_starts=col_starts)
        elif isinstance(r, mfem.par.HypreParMatrix):
            self[0] = r
        elif r is None:
            self[0] = r
        else:
            assert False, "unknonw matrix element"
        if isinstance(i, csr_matrix):
            self[1] = ToHypreParCSR(i, col_starts=col_starts)
        elif isinstance(i, mfem.par.HypreParMatrix):
            self[1] = i
        elif i is None:
            self[1] = i
        else:
            assert False, "unknonw matrix element"

    def GetOwns(self):
        flags = [None] * 6
        if self[0] is not None:
            flags[0] = self[0].OwnsDiag()
            flags[1] = self[0].OwnsOffd()
            flags[2] = self[0].OwnsColMap()
        if self[1] is not None:
            flags[3] = self[1].OwnsDiag()
            flags[4] = self[1].OwnsOffd()
            flags[5] = self[1].OwnsColMap()
        return flags

    def __repr__(self):
        return "CHypreMat" + str(self.shape) + "[" + str(self.lshape) + "]"

    def isComplex(self):
        return not (self[1] is None)

    def __mul__(self, other):  # A * B or A * v
        if isinstance(other, CHypreMat):
            return CHypreMat(*ParMultComplex(self, other))
        elif isinstance(other, CHypreVec):
            v = CHypreVec(*ParMultVecComplex(self, other))
            v._horizontal = other._horizontal
            return v
        elif isinstance(other, Number):
            if np.iscomplexobj(other):
                other = complex(other)
                r = other.real
                i = other.imag
            else:
                r = other
                i = 0
            if self[0] is not None and self[1] is not None:
                R = mfem.par.Add(r, self[0], -i, self[1])
                I = mfem.par.Add(i, self[0], r, self[1])
            elif self[0] is not None:
                #R = mfem.par.Add(r, self[0], 0.0, self[0])
                R = mfem.par.Add(r, self[0], 0.0, self[0])
                if i != 0.0:
                    I = mfem.par.Add(i, self[0], 0.0, self[0])
                else:
                    I = None
            elif self[1] is not None:
                if r != 0.0:
                    I = mfem.par.Add(r, self[1], 0.0, self[1])
                else:
                    I = None
                R = mfem.par.Add(-i, self[1], 0.0, self[1])
            else:
                assert False, "this mode is not supported"
                R = None
                I = None
            if R is not None:
                R.CopyRowStarts()
                R.CopyColStarts()
            if I is not None:
                I.CopyRowStarts()
                I.CopyColStarts()
            return CHypreMat(R, I)
        raise ValueError("argument should be CHypreMat/Vec")

    def __rmul__(self, other):
        if not isinstance(other, CHypreMat):
            raise ValueError(
                "argument should be CHypreMat")
        return CHypreMat(*ParMultComplex(other, self))

    def __add__(self, other):  # A + B
        if not isinstance(other, CHypreMat):
            raise ValueError(
                "argument should be CHypreMat")
        if self[0] is not None and other[0] is not None:
            r = ParAdd(self[0], other[0])
        elif self[0] is not None:
            r = ToHypreParCSR(ToScipyCoo(self[0]).tocsr())
        elif other[0] is not None:
            r = ToHypreParCSR(ToScipyCoo(other[0]).tocsr())
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            i = ParAdd(self[1], other[1])
#            i =  mfem.par.add_hypre(1.0, self[1], 1.0, other[1])
        elif self[1] is not None:
            i = ToHypreParCSR(ToScipyCoo(self[1]).tocsr())
        elif other[1] is not None:
            i = ToHypreParCSR(ToScipyCoo(other[1]).tocsr())
        else:
            i = None

        if r is not None:
            r.CopyRowStarts()
            r.CopyColStarts()
        if i is not None:
            i.CopyRowStarts()
            i.CopyColStarts()

        return CHypreMat(r, i)

    def __sub__(self, other):  # A - B
        if not isinstance(other, CHypreMat):
            raise ValueError(
                "argument should be CHypreMat")
        if self[0] is not None and other[0] is not None:
            other[0] *= -1
            r = ParAdd(self[0], other[0])
            other[0] *= -1
        elif self[0] is not None:
            r = ToHypreParCSR(ToScipyCoo(self[0]).tocsr())
        elif other[0] is not None:
            r = mfem.par.Add(-1, othe[0], 0.0, other[0])
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            other[1] *= -1
            i = ParAdd(self[1], other[1])
            other[1] *= -1
        elif self[1] is not None:
            i = ToHypreParCSR(ToScipyCoo(self[1]).tocsr())
        elif other[1] is not None:
            i = mfem.par.Add(-1, othe[1], 0.0, other[1])
        else:
            i = None

        if r is not None:
            r.CopyRowStarts()
            r.CopyColStarts()
        if i is not None:
            i.CopyRowStarts()
            i.CopyColStarts()

        return CHypreMat(r, i)

    def __neg__(self):  # -B
        r = None
        i = None
        if self[0] is not None:
            r = mfem.par.Add(-1, self[0], 0.0, self[0])
            r.CopyRowStarts()
            r.CopyColStarts()
        if self[1] is not None:
            i = mfem.par.Add(-1, self[1], 0.0, self[0])
            i.CopyRowStarts()
            i.CopyColStarts()

        return CHypreMat(r, i)

    @property
    def imag(self):
        return self[1]

    @imag.setter
    def imag(self, value):
        self[1] = value

    @property
    def real(self):
        return self[0]

    @real.setter
    def real(self, value):
        self[0] = value

    def GetColPartArray(self):
        if self[0] is not None:
            return self[0].GetColPartArray()
        else:
            return self[1].GetColPartArray()

    def GetRowPartArray(self):
        if self[0] is not None:
            return self[0].GetRowPartArray()
        else:
            return self[1].GetRowPartArray()

    def __imul__(self, other):
        if self[0] is not None:
            self[0] *= other
        if self[1] is not None:
            self[1] *= other
        return self

    def dot(self, other):
        return self.__mul__(other)

    def transpose(self):
        '''
        transpose of matrix
        this method returns a new matrix
        '''
        R = self[0].Transpose() if self[0] is not None else None
        I = self[1].Transpose() if self[1] is not None else None
        return CHypreMat(R, I)
        # return CHypreMat(self[0].Transpose(), self[1].Transpose())

    def conj(self, inplace=True):
        '''
        complex conjugate
          if copy is on, imaginary part becomes different object
        '''
        if self[1] is None:
            return self

        if not inplace:
            self[1] = ToHypreParCSR(-ToScipyCoo(self[1]).tocsr())
        else:
            self[1] *= -1.0
        return self

    def rap(self, B):
        '''
        B^* A B
        '''
        ret = B.conj().transpose() * self
        ret = ret * (B.conj())
        # this should put B back to original
        return ret

    def setDiag(self, idx, value=1.0):
        if self[0] is not None:
            self[0] = ResetHypreDiag(self[0], idx, value=np.real(value))
        if self[1] is not None:
            self[1] = ResetHypreDiag(self[1], idx, value=np.imag(value))

    def getDiag(self, idx):
        if self[0] is not None:
            diagvalue = ReadHypreDiag(self[0], idx)
        else:
            diagvalue = complex(0.0)

        if self[1] is not None:
            diagvalue = diagvalue + 1j*ReadHypreDiag(self[1], idx)
        else:
            diagvalue = diagvalue + 0j
        return diagvalue

    def resetCol(self, idx, inplace=True):
        if self[0] is not None:
            r = ResetHypreCol(self[0], idx)
        else:
            r = None
        if self[1] is not None:
            i = ResetHypreCol(self[1], idx)
        else:
            i = None
        if inplace:
            self[0] = r
            self[1] = i
            return self
        else:
            return CHypreMat(r, i)

    def resetRow(self, idx, inplace=True):
        if self[0] is not None:
            r = ResetHypreRow(self[0], idx)
        else:
            r = None
        if self[1] is not None:
            i = ResetHypreRow(self[1], idx)
        else:
            i = None
        if inplace:
            self[0] = r
            self[1] = i
            return self
        else:
            return CHypreMat(r, i)

    def selectRows(self, idx):
        '''
        idx is global index
        '''
        rpart = self[0].GetRowPartArray()
        #rpart[2] = self[0].GetGlobalNumRows()
        cpart = self[0].GetColPartArray()
        #cpart[2] = self[0].GetGlobalNumCols()

        idx = idx[idx >= rpart[0]]
        idx = idx[idx < rpart[1]]
        idx = idx - rpart[0]

        csr = ToScipyCoo(self[0]).tocsr()

        csr = csr[idx, :]
        r = ToHypreParCSR(csr, col_starts=cpart)
        if self[1] is not None:
            csr = ToScipyCoo(self[1]).tocsr()
            csr = csr[idx, :]
            i = ToHypreParCSR(csr, col_starts=cpart)
        else:
            i = None
        return CHypreMat(r, i)

    def selectCols(self, idx):
        '''
        idx is global index
        '''
        cpart = self[0].GetColPartArray()
        #cpart[2] = self[0].GetGlobalNumCols()
        rpart = self[0].GetRowPartArray()
        #rpart[2] = self[0].GetGlobalNumRows()

        idx = idx[idx >= cpart[0]]
        idx = idx[idx < cpart[1]]
        idx = idx - cpart[0]

        mat = self.transpose()
        csr = ToScipyCoo(mat[0]).tocsr()
        csr = csr[idx, :]
        r = ToHypreParCSR(csr, col_starts=rpart)
        if self[1] is not None:
            csr = ToScipyCoo(mat[1]).tocsr()
            csr = csr[idx, :]
            i = ToHypreParCSR(csr, col_starts=rpart)
        else:
            i = None
        mat = CHypreMat(r, i).transpose()
        '''
        if (cpart == rpart).all():
           csr = ToScipyCoo(mat[0]).tocsr()
           mat[0] = ToHypreParCSR(csr, col_starts =rpart)
           if mat[1] is not None:
              csr = ToScipyCoo(mat[1]).tocsr()
              mat[1] = ToHypreParCSR(csr, col_starts =rpart)
        '''
        return mat

    @property
    def nnz(self):
        if self[0] is not None and self[1] is not None:
            return self[0].NNZ(), self[1].NNZ()
        if self[0] is not None:
            return self[0].NNZ()
        if self[1] is not None:
            return self[1].NNZ()

    def true_nnz(self):
        '''
        more expensive version which reports nnz after
        eliminating all zero entries
        '''
        #coo = self.get_local_coo()
        if self[0] is not None:
            nnz0, tnnz0 = self[0].get_local_true_nnz()
        else:
            nnz0 = 0
            tnnz0 = 0
        if self[1] is not None:
            nnz1, tnnz1 = self[1].get_local_true_nnz()
        else:
            nnz1 = 0
            tnnz1 = 0
        # print nnz0, tnnz0, nnz1, tnnz1
        return tnnz0, tnnz1

    def m(self):
        '''
        return global row number: two number should be the same
        '''
        ans = []
        if self[0] is not None:
            ans.append(self[0].M())
        if self[1] is not None:
            ans.append(self[1].M())

        if len(ans) == 2 and ans[0] != ans[1]:
            raise ValueError(
                "data format error, real and imag should have same size")
        return ans[0]

    def n(self):
        '''
        return global col number: two number should be the same
        '''
        ans = []
        if self[0] is not None:
            ans.append(self[0].N())
        if self[1] is not None:
            ans.append(self[1].N())
        if len(ans) == 2 and ans[0] != ans[1]:
            raise ValueError(
                "data format error, real and imag should have same size")
        return ans[0]

    @property
    def shape(self):
        if self[0] is not None:
            return (self[0].GetGlobalNumRows(), self[0].GetGlobalNumCols())
        elif self[1] is not None:
            return (self[1].GetGlobalNumRows(), self[1].GetGlobalNumCols())
        else:
            return (0, 0)

    @property
    def lshape(self):
        if self[0] is not None:
            return (self[0].GetNumRows(), self[0].GetNumCols())
        elif self[1] is not None:
            return (self[1].GetNumRows(), self[1].GetNumCols())
        else:
            return (0, 0)

    def get_local_coo(self):
        if self.isComplex():
            return (ToScipyCoo(self[0]) + 1j * ToScipyCoo(self[1])).tocoo()
        else:
            return ToScipyCoo(self[0])

    def get_global_coo(self):
        lcoo = self.get_local_coo()
        gcoo = coo_matrix(self.shape)
        gcoo.data = lcoo.data

        gcoo.row = lcoo.row + self.GetRowPartArray()[0]
        gcoo.col = lcoo.col
        return gcoo

    def get_squaremat_from_right(self):
        '''
        squre matrix which can be multipled from the right of self.
        '''
        size = self.shape[1]
        if self[0] is not None:
            part = self[0].GetColPartArray()
            width = self[0].GetGlobalNumCols()
        elif self[1] is not None:
            part = self[1].GetColPartArray()
            width = self[1].GetGlobalNumCols()
        else:
            raise ValueError("CHypreMat is empty")
        return SquareCHypreMat(width, part, real=(self[1] is None))

    def elimination_matrix(self, idx):
        #  # local version
        #  ret = lil_matrix((len(nonzeros), self.shape[0]))
        #  for k, z in enumerate(nonzeros):
        #     ret[k, z] = 1.
        #  return ret.tocoo()
        cpart = self.GetColPartArray()
        rpart = self.GetRowPartArray()

        idx = idx[idx >= rpart[0]]
        idx = idx[idx < rpart[1]]
        idx = idx - rpart[0]

        shape = (len(idx), self.shape[1])
        # print shape, idx + rpart[0]
        elil = lil_matrix(shape)
        for i, j in enumerate(idx):
            elil[i, j + rpart[0]] = 1

        r = ToHypreParCSR(elil.tocsr(), col_starts=cpart)
        return CHypreMat(r, None)

    def eliminate_RowsCols(self, B, tdof, inplace=False, diagpolicy=0):
        # note: tdof is a valued viewed in each node. MyTDoF offset is
        # subtracted
        tdof1 = mfem.par.intArray(tdof)
        if not inplace:
            #print("inplace flag off copying ParCSR")
            if self[0] is not None:
                r = ToHypreParCSR(ToScipyCoo(self[0]).tocsr())

            else:
                r = None
            if self[1] is not None:
                i = ToHypreParCSR(ToScipyCoo(self[1]).tocsr())

            else:
                i = None
            target = CHypreMat(r, i)
        else:
            target = self

        row0 = self.GetRowPartArray()[0]
        gtdof = list(np.array(tdof, dtype=int) + row0)

        if diagpolicy == 0:
            diagAe = target.getDiag(gtdof) - 1
            diagA = 1.0
        else:
            diagAe = 0.0
            diagA = target.getDiag(gtdof)

        if target[0] is not None:
            Aer = target[0].EliminateRowsCols(tdof1)
            Aer.CopyRowStarts()
            Aer.CopyColStarts()
            #row0 = Aer.GetRowPartArray()[0]
        else:
            Aer = None

        if target[1] is not None:
            Aei = target[1].EliminateRowsCols(tdof1)
            Aei.CopyRowStarts()
            Aei.CopyColStarts()
            #row0 = Aei.GetRowPartArray()[0]
        else:
            Aei = None

        Ae = CHypreMat(Aer, Aei)

        target.setDiag(gtdof, value=diagA)
        Ae.setDiag(gtdof, value=diagAe)

        # if diagpolicy == 0:
        #    part = B.GetPartitioningArray()
        #    xxx = np.ones(part[1] - part[0], dtype=complex)
        #    xxx[tdof] = diagA
        #    B *= xxx
        # if diagpolicy == 1:
        #    print(tdof)
        #    print(diagA.shape)
        #    diagA = diagA[tdof]

        B.set_elements(gtdof, diagA)

        return Ae, target, B

    @property
    def isHypre(self):
        return True


def SquareCHypreMat(width, part, real=False):
    from scipy.sparse import csr_matrix
    lr = part[1] - part[0]
    m1 = csr_matrix((lr, width))
    if real:
        m2 = None
    else:
        m2 = csr_matrix((lr, width))
    return CHypreMat(m1, m2)


def Array2CHypreVec(array, part=None, horizontal=False):
    '''
    convert array in rank =0 to
    distributed Hypre 1D Matrix (size = m x 1)
    '''
    from mpi4py import MPI
    isComplex = MPI.COMM_WORLD.bcast(np.iscomplexobj(array), root=0)
    if isComplex:
        if array is None:
            rarray = None
            iarray = None
        else:
            rarray = array.real
            iarray = array.imag
        return CHypreVec(Array2HypreVec(rarray, part),
                         Array2HypreVec(iarray, part), horizontal=horizontal)
    else:
        if array is None:
            rarray = None
        else:
            rarray = array
        return CHypreVec(Array2HypreVec(rarray, part),
                         None, horizontal=horizontal)


def CHypreVec2Array(array):
    from mpi4py import MPI
    myid = MPI.COMM_WORLD.rank

    if array[0] is not None:
        r = HypreVec2Array(array[0])
    else:
        if myid == 0:
            r = 0.0
        else:
            r = None
    if array[1] is None:
        return r
    else:
        i = HypreVec2Array(array[1])
    if i is None:
        return r
    else:
        if myid == 0:
            return r + 1j * i
        else:
            return None


def CHypreMat2Coo(mat):
    print("CHYPREMat2Coo: deprecated,  Use class method !!!!")
    if mat.isComplex():
        return ToScipyCoo(mat.real) + 1j * ToScipyCoo(mat.imag)
    else:
        return ToScipyCoo(mat.real)


def LF2PyVec(rlf, ilf=None, horizontal=False):
    if MFEM_PAR:
        '''
        From ParLF to CHypreVec
        '''
        rv = rlf.ParallelAssemble()
        rv.thisown = True
        if ilf is not None:
            iv = ilf.ParallelAssemble()
            iv.thisown = True
        else:
            iv = None
        return CHypreVec(rv, iv, horizontal=horizontal)
    else:
        b1 = rlf.GetDataArray().copy()  # ; rlf.thisown = False
        if ilf is not None:
            b2 = ilf.GetDataArray()  # ; ilf.thisown = False
            b1 = b1 + 1j * b2
        if horizontal:
            return b1.reshape((1, -1))
        else:
            return b1.reshape((-1, 1))


LinearForm2PyVector = LF2PyVec


def MfemVec2PyVec(rlf, ilf=None, horizontal=False):
    b1 = rlf.GetDataArray().copy()  # ; rlf.thisown = False
    if ilf is not None:
        b2 = ilf.GetDataArray()  # ; ilf.thisown = False
    else:
        b2 = None

    if MFEM_PAR:
        b1 = ToHypreParVec(b1)
        if b2 is not None:
            b2 = ToHypreParVec(b2)
        return CHypreVec(b1, b2, horizontal=horizontal)
    else:
        if b2 is not None:
            b1 = b1 + 1j * b2
        if horizontal:
            return b1.reshape((1, -1))
        else:
            return b1.reshape((-1, 1))


def Array2PyVec(array, part=None, horizontal=False):
    '''
    convert array in rank = 0  to
    distributed Hypre 1D Matrix (size = m x 1)
    '''

    if MFEM_PAR:
        return Array2CHypreVec(array, part=part, horizontal=horizontal)
    else:
        if horizontal:
            return array.reshape((1, -1))
        else:
            return array.reshape((-1, 1))


def BF2PyMat(rbf, ibf=None, finalize=False):
    '''
    Convert pair of BilinearForms to CHypreMat or
    ScipySparsematrix
    '''
    if finalize:
        rbf.Finalize()
        if ibf is not None:
            ibf.Finalize()

    if MFEM_PAR:
        M1 = rbf.ParallelAssemble()
        M1.thisown = True
        if ibf is not None:
            M2 = ibf.ParallelAssemble()
            M2.thisown = True
        else:
            M2 = None
        return CHypreMat(M1, M2)
    else:
        from mfem.common.sparse_utils import sparsemat_to_scipycsr
        M1 = rbf.SpMat()
        if ibf is None:
            return sparsemat_to_scipycsr(M1, dtype=float)
        if ibf is not None:
            M2 = ibf.SpMat()
            m1 = sparsemat_to_scipycsr(M1, dtype=float).tolil()
            m2 = sparsemat_to_scipycsr(M2, dtype=complex).tolil()
            m = m1 + 1j * m2
            m = m.tocsr()
        return m


BilinearForm2PyMatix = BF2PyMat


def MfemMat2PyMat(M1, M2=None):
    '''
    Convert pair of SpMat/HypreParCSR to CHypreMat or
    ScipySparsematrix. This is simpler version of BF2PyMat, only
    difference is it skippes convertion from BF to Matrix.
    '''
    from mfem.common.sparse_utils import sparsemat_to_scipycsr
    if MFEM_PAR:
        return CHypreMat(M1, M2)
    else:
        if M2 is None:
            return sparsemat_to_scipycsr(M1, dtype=float)
        else:
            m1 = sparsemat_to_scipycsr(M1, dtype=float).tolil()
            m2 = sparsemat_to_scipycsr(M2, dtype=complex).tolil()
            m = m1 + 1j * m2
            m = m.tocsr()
            return m


def EmptySquarePyMat(m, col_starts=None):
    from scipy.sparse import csr_matrix
    if MFEM_PAR:
        if col_starts is None:
            col_starts = get_assumed_patitioning(m)
        rows = col_starts[1] - col_starts[0]
        m1 = csr_matrix((rows, m))
        return CHypreMat(m1, None, )
    else:
        from scipy.sparse import csr_matrix
        return csr_matrix((m, m))


def IdentityPyMat(m, col_starts=None, diag=1.0):
    from scipy.sparse import coo_matrix, lil_matrix
    if MFEM_PAR:
        if col_starts is None:
            col_starts = get_assumed_patitioning(m)
        rows = col_starts[1] - col_starts[0]

        if np.iscomplexobj(diag):
            real = diag.real
            imag = diag.imag
        else:
            real = float(np.real(diag))
            imag = 0.0
        # if real != 0.0:
        m1 = lil_matrix((rows, m))
        for i in range(rows):
            m1[i, i + col_starts[0]] = real
        m1 = m1.tocsr()

        if imag != 0.0:
            m2 = lil_matrix((rows, m))
            for i in range(rows):
                m2[i, i + col_starts[0]] = imag
            m2 = m2.tocsr()
        else:
            m2 = None
        return CHypreMat(m1, m2, )
    else:
        m1 = coo_matrix((m, m))
        m1.setdiag(np.zeros(m) + diag)
        return m1.tocsr()


def HStackPyVec(vecs, col_starts=None):
    '''
    horizontally stack vertical vectors to generate
    PyMat
    '''
    from scipy.sparse import csr_matrix
    if MFEM_PAR:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rows = vecs[0].GetPartitioningArray()
        if col_starts is None:
            col_starts = get_assumed_patitioning(len(vecs))

        isComplex = any([v.isComplex() for v in vecs])
        mat = np.hstack([np.atleast_2d(v.toarray()).transpose() for v in vecs])
        if isComplex:
            m1 = csr_matrix(mat.real)
            m2 = csr_matrix(mat.imag)
        else:
            m1 = csr_matrix(mat)
            m2 = None
        ret = CHypreMat(m1, m2, col_starts=col_starts)
        return ret
    else:
        return csr_matrix(np.hstack(vecs))


def PyVec2PyMat(vec, col_starts=None):
    from scipy.sparse import csr_matrix
    if MFEM_PAR:
        '''
        vec must be vertical
        '''
        assert not vec._horizontal, "PyVec must be vertical"

        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rows = vec.GetPartitioningArray()
        if col_starts is None:
            col_starts = get_assumed_patitioning(1)

        isComplex = vec.isComplex()
        mat = np.atleast_2d(vec.toarray()).transpose()

        if isComplex:
            m1 = csr_matrix(mat.real)
            m2 = csr_matrix(mat.imag)
        else:
            m1 = csr_matrix(mat)
            m2 = None
        # print m1.shape
        ret = CHypreMat(m1, m2, col_starts=col_starts)
        # print "returning", ret,
        return ret
    else:
        return csr_matrix(vec)
