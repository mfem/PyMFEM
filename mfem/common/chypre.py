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
from scipy.sparse  import csr_matrix, coo_matrix, lil_matrix
from mfem.common.parcsr_extra import *

# DO NOT IMPORT MPI in Global, sicne some routins will be used
# in serial mode too.
try:
   import mfem.par
   MFEM_PAR = True
except:
   MFEM_PAR = False

class CHypreVec(list):
    def __init__(self, r = None, i = None, horizontal = False):
        list.__init__(self, [None]*2)
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
        return "CHypreVec"+str(self.shape)+"["+str(part[1]-part[0])+"]"
            
    @property
    def imag(self):
        return self[1]
    @imag.setter
    def imag(self, value):
        self[1]  = value
    @property
    def real(self):
        return self[0]
    @real.setter
    def real(self, value):
        self[0]  = value
    @property
    def shape(self):
        if self._horizontal:
            return 1, self[0].GlobalSize()
        else:
            return self[0].GlobalSize(), 1
            
    def isComplex(self):
        return not (self[1] is None)
     
    def GetPartitioningArray(self):
        if self[0] is not None:
            part = self[0].GetPartitioningArray()
            part[2] = self[0].GlobalSize()
        elif self[1] is not None:         
            prat =  self[1].GetPartitioningArray()
            part[2] = self[1].GlobalSize()         
        else:
            raise ValueError("CHypreVec is empty")
        return part
     
    def __imul__(self, other):
        if self[0] is not None:
            self[0] *= other
        if self[1] is not None:
            self[1] *= other
        return self

    def __add__(self, other):
        Vector = mfem.par.HypreParVector
        if self[0] is not None and other[0] is not None:
            r = Vector(self[0])            
            add_vector(self[0], others[0], r)
        elif self[0] is not None:
            r = Vector(self[0])
        elif other[0] is not None:
            r = Vector(other[0])
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            i = Vector(self[1])            
            add_vector(self[1], others[1], i)
        elif self[1] is not None:
            i = Vector(self[1])
        elif other[1] is not None:
            i = Vector(other[1])
        else:
            i = None
        return CHypreVec(r, i, horizontal = self._horizontal)

    def __sub__(self, other):
        from mfem.par import HypreParVector as Vector
        from mfem.par import add_vector
        add_vector 
        if self[0] is not None and other[0] is not None:
            r = Vector(self[0])            
            add_vector(self[0], -1, other[0], r)
        elif self[0] is not None:
            r = Vector(self[0])
        elif other[0] is not None:
            r = Vector(other[0])
            r *=-1.
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            i = Vector(self[1])            
            add_vector(self[1], -1, other[1], i)
        elif self[1] is not None:
            i = Vector(self[1])
        elif other[1] is not None:
            i = Vector(other[1])
            i *=-1.
        else:
            i = None

        return CHypreVec(r, i, horizontal = self._horizontal)
    
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

    
    def set_element(self, i, v):
        part = self.GetPartitioningArray()
        if part[0] <= i and i < part[1]:
            v = complex(v)
            if self[0] is not None:
                self[0][int(i - part[0])] = v.real
            if self[1] is not None:            
                self[1][int(i - part[0])] = v.imag
                
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
            
    def get_squaremat_from_right(self):
        '''
        squre matrix which can be multipled from the right of self.
        '''
        if not self._horizontal:
            raise ValueError("Vector orientation is not right")
        
        part = self.GetPartitioningArray()
        return SquareCHypreMat(part)
    
    def transpose(self):
        self._horizontal = not self.horizontal

    def _do_reset(self, v, idx):
        # ownership is transferrd to a new vector.
        ownership = v.GetOwnership()
        data = v.GetDataArray()
        part = v.GetPartitioningArray()
        for i in idx:
            if i >= part[1]:continue                
            if i < part[0]:continue
            data[i - part[0]] = 0
            
        ret = ToHypreParVec(data)
        ret.SetOwnership(ownership)
        v.SetOwnership(0)
        return  ret
        
    def resetCol(self, idx):
        if self._horizontal:
            if self[0] is not None:
                self[0] = self._do_reset(self[0], idx)
            if self[1] is not None:
                self[1] = self._do_reset(self[1], idx)
        else:            
            if 0 in idx: self *= 0.0
            
    def resetRow(self, idx):
        if self._horizontal:
            if 0 in idx: self *= 0.0
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
        return  ret

    def selectRows(self, idx):
        '''
        idx is global index
        '''
        if self._horizontal:
            if not 0 in idx: raise ValueError("VectorSize becomes zero")
            return self

        part = self.GetPartitioningArray()                        
        idx = idx[idx >= part[0]]
        idx = idx[idx < part[1]]
        lidx = idx - part[0]

        r = None; i = None        
        if not self._horizontal:
            if self[0] is not None:
                r = self._do_select(self[0], lidx)
            if self[1] is not None:
                i = self._do_select(self[1], lidx)
        return CHypreVec(r, i, horizontal = self._horizontal)
     
    def selectCols(self, idx):
        '''
        idx is global index
        '''
        if not self._horizontal:
            if not 0 in idx: raise ValueError("VectorSize becomes zero")
            return self
         
        part = self.GetPartitioningArray()                        
        idx = idx[idx >= part[0]]
        idx = idx[idx < part[1]]
        lidx = idx - part[0]

        r = None; i = None
        if self._horizontal:
            if self[0] is not None:
                r = self._do_select(self[0], lidx)
            if self[1] is not None:
                i = self._do_select(self[1], lidx)
        return CHypreVec(r, i, horizontal = self._horizontal)
     
    def GlobalVector(self):
        if self[0] is not None:
            gv = self[0].GlobalVector()
            v = gv.GetDataArray().copy()  
        else:
            v = 0.0
        if self[1] is not None:            
            gv = self[1].GlobalVector()
            v = v + 1j*gv.GetDataArray()
        return v

    def toarray(self):
        '''
        numpy array of local vector
        '''
        if self[0] is not None:
            v = self[0].GetDataArray()
        else:
            v = 0.0
        if self[1] is not None:            
            v = v + 1j*self[1].GetDataArray()
        return v

    def isAllZero(self):
        return any(self.GlobalVector())
     
    def get_global_coo(self):
        data = self[0].GetDataArray()
        if self[1] is not None:
            data = data + 1j*self[1].GetDataArray()
        gcoo = coo_matrix(self.shape, dtype = data.dtype)
        gcoo.data = data
        part = self.GetPartitioningArray()        
        if self._horizontal:
            gcoo.row  = np.zeros(len(data))
            gcoo.col  = np.arange(len(data))+ part[0]
        else:
            gcoo.col  = np.zeros(len(data))
            gcoo.row  = np.arange(len(data))+ part[0]
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
            data = data + 1j*self[1].GetDataArray()
        local_nnz = np.sum(data == 0.)

        from mpi4py import MPI
        comm     = MPI.COMM_WORLD     
        nnz = comm.allgather(local_nnz)
        return nnz
     
    @property     
    def isHypre(self):
        return True

class CHypreMat(list):
    def __init__(self, r = None, i = None, col_starts = None):
        list.__init__(self, [None]*2)
        if isinstance(r, csr_matrix):
            self[0] = ToHypreParCSR(r, col_starts = col_starts)
        else:
            self[0] = r
        if isinstance(i, csr_matrix):
            self[1] = ToHypreParCSR(i, col_starts = col_starts)
        else:
            self[1] = i

    def __repr__(self):
        return "CHypreMat"+str(self.shape) +"["+str(self.lshape)+"]"
    
    def isComplex(self):
        return not (self[1] is None)

    def __mul__(self, other): # A * B or A * v
        if isinstance(other, CHypreMat):
            return CHypreMat(*ParMultComplex(self, other))
        if isinstance(other, CHypreVec):
            v =  CHypreVec(*ParMultVecComplex(self, other))
            v._horizontal = other._horizontal
            return v
        raise ValueError(
                   "argument should be CHypreMat/Vec")
    
    def __rmul__(self, other):
        if not isinstance(other, CHypreMat):
             raise ValueError(
                   "argument should be CHypreMat")
        return CHypreMat(*ParMultComplex(other, self))
        
    def __add__(self, other): #A + B
        if not isinstance(other, CHypreMat):
             raise ValueError(
                   "argument should be CHypreMat")
        if self[0] is not None and other[0] is not None:            
            r =  ParAdd(self[0], other[0])
        elif self[0] is not None:
            r = self[0]
        elif other[0] is not None:
            r = other[0]
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            i =  ParAdd(self[1], other[1])           
#            i =  mfem.par.add_hypre(1.0, self[1], 1.0, other[1])
        elif self[0] is not None:
            i = self[1]
        elif other[0] is not None:
            i = other[1]
        else:
            i = None

        return CHypreMat(r, i)        

    def __sub__(self, other): #A - B
        if not isinstance(other, CHypreMat):
             raise ValueError(
                   "argument should be CHypreMat")
        if self[0] is not None and other[0] is not None:
            other[0] *= -1
            r =  ParAdd(self[0], other[0])
            other[0] *= -1            
        elif self[0] is not None:
            r = self[0]
        elif other[0] is not None:
            r = other[0]
        else:
            r = None
        if self[1] is not None and other[1] is not None:
            other[1] *= -1           
            i =  ParAdd(self[1], other[1])
            other[1] *= -1                       
        elif self[0] is not None:
            i = self[1]
        elif other[0] is not None:
            i = other[1]
        else:
            i = None

        return CHypreMat(r, i)        
       
    def __neg__(self): #-B
       r = None; i = None
       if self[0] is not None:
           r = ToHypreParCSR((- ToScipyCoo(self[0])).tocsr())
       if self[1] is not None:
           i = ToHypreParCSR((- ToScipyCoo(self[1])).tocsr())
       return CHypreMat(r, i)


    @property
    def imag(self):
        return self[1]
    @imag.setter
    def imag(self, value):
        self[1]  = value
    @property
    def real(self):
        return self[0]
    @real.setter
    def real(self, value):
        self[0]  = value

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
        #return CHypreMat(self[0].Transpose(), self[1].Transpose())

    def conj(self, inplace = True):
        '''
        complex conjugate
          if copy is on, imaginary part becomes different object
        '''
        if self[1] is None: return self
        
        if not inplace:
           self[1] = ToHypreParCSR(-ToScipyCoo(self[1]).tocsr())
        else:
           self[1] *= -1.0
        return self

    def rap(self, B):
        '''
        B^* A B
        '''
        ret = B.conj().transpose()*self
        ret = ret * (B.conj())
        # this should put B back to original
        return ret

    def setDiag(self, idx, value = 1.0):
        if self[0] is not None:
            self[0] = ResetHypreDiag(self[0], idx, value = float(value))
        if self[1] is not None:            
            self[1] = ResetHypreDiag(self[1], idx, value = float(np.imag(value)))
            
    def resetCol(self, idx):
        if self[0] is not None:
            self[0] = ResetHypreCol(self[0], idx)
        if self[1] is not None:            
            self[1] = ResetHypreCol(self[1], idx)
            
    def resetRow(self, idx):
        if self[0] is not None:
            self[0] = ResetHypreRow(self[0], idx)
        if self[1] is not None:            
            self[1] = ResetHypreRow(self[1], idx)
            
    def selectRows(self, idx):
        '''
        idx is global index
        '''
        rpart = self[0].GetRowPartArray()
        rpart[2] = self[0].GetGlobalNumRows()            
        cpart = self[0].GetColPartArray()
        cpart[2] = self[0].GetGlobalNumCols()
        
        idx = idx[idx >= rpart[0]]
        idx = idx[idx < rpart[1]]
        idx = idx - rpart[0]

        csr = ToScipyCoo(self[0]).tocsr()

        csr = csr[idx, :]
        r = ToHypreParCSR(csr, col_starts = cpart)
        if self[1] is not None:
           csr = ToScipyCoo(self[1]).tocsr()
           csr = csr[idx, :]
           i = ToHypreParCSR(csr, col_starts =cpart)
        else:
           i = None
        return CHypreMat(r, i)
       
    def selectCols(self, idx):
        '''
        idx is global index
        '''
        cpart = self[0].GetColPartArray()
        cpart[2] = self[0].GetGlobalNumCols()
        rpart = self[0].GetRowPartArray()
        rpart[2] = self[0].GetGlobalNumRows()
        
        idx = idx[idx >= cpart[0]]
        idx = idx[idx < cpart[1]]
        idx = idx - cpart[0]

        mat = self.transpose()        
        csr = ToScipyCoo(mat[0]).tocsr()
        csr = csr[idx, :]
        r = ToHypreParCSR(csr, col_starts = rpart)
        if self[1] is not None:
           csr = ToScipyCoo(mat[1]).tocsr()
           csr = csr[idx, :]
           i = ToHypreParCSR(csr, col_starts =rpart)
        else:
           i = None
        mat = CHypreMat(r, i).transpose()
        return mat
    
    @property
    def nnz(self):
        if self[0] is not None and self[1] is not None:
            return self[0].NNZ(), self[1].NNZ()
        if self[0] is not None: return self[0].NNZ()
        if self[1] is not None: return self[1].NNZ()

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
        #print nnz0, tnnz0, nnz1, tnnz1
        return tnnz0, tnnz1

    def m(self):
        '''
        return global row number: two number should be the same
        '''
        ans = []
        if self[0] is not None: ans.append(self[0].M())
        if self[1] is not None: ans.append(self[1].M())

        if len(ans) == 2 and ans[0] != ans[1]:
            raise ValueError("data format error, real and imag should have same size")
        return ans[0]
    
    def n(self):
        '''
        return global col number: two number should be the same
        '''
        ans = []
        if self[0] is not None: ans.append(self[0].N())
        if self[1] is not None: ans.append(self[1].N())
        if len(ans) == 2 and ans[0] != ans[1]:
            raise ValueError("data format error, real and imag should have same size")
        return ans[0]
    @property
    def shape(self):
        if self[0] is not None:
            return (self[0].GetGlobalNumRows(), self[0].GetGlobalNumCols())
        elif self[1] is not None:        
            return (self[1].GetGlobalNumRows(), self[1].GetGlobalNumCols())
        else:
            return (0,0)
    @property
    def lshape(self):
        if self[0] is not None:
            return (self[0].GetNumRows(), self[0].GetNumCols())
        elif self[1] is not None:        
            return (self[1].GetNumRows(), self[1].GetNumCols())
        else:
            return (0,0)

    def get_local_coo(self):
        if self.isComplex():
            return (ToScipyCoo(self[0]) + 1j*ToScipyCoo(self[1])).tocoo()
        else:
            return ToScipyCoo(self[0])

    def get_global_coo(self):
        lcoo = self.get_local_coo()
        gcoo = coo_matrix(self.shape)
        gcoo.data = lcoo.data

        gcoo.row  = lcoo.row + self.GetRowPartArray()[0]
        gcoo.col  = lcoo.col
        return gcoo
        
    def get_squaremat_from_right(self):
        '''
        squre matrix which can be multipled from the right of self.
        '''
        size = self.shape[1]
        if self[0] is not None:
            part = self[0].GetColPartArray()
            part[2] = self[0].GetGlobalNumCols()
        elif self[1] is not None:
            part = self[1].GetColPartArray()
            part[2] = self[1].GetGlobalNumCols()            
        else:
            raise ValueError("CHypreMat is empty")
        return SquareCHypreMat(part)
     
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
        print shape, idx + rpart[0]
        elil =  lil_matrix(shape)
        for i, j in enumerate(idx):
           elil[i, j + rpart[0]] = 1
       
        r = ToHypreParCSR(elil.tocsr(), col_starts =cpart)
        return CHypreMat(r, None)
     
    @property     
    def isHypre(self):
        return True
    
def SquareCHypreMat(part):
    from scipy.sparse import csr_matrix
    lr = part[1]-part[0]
    c  = part[2]
    m1 = csr_matrix((lr, c))
    m2 = csr_matrix((lr, c))
    return CHypreMat(m1, m2)

def Array2CHypreVec(array, part): 
    '''
    convert array in rank (default = 0)  to 
    distributed Hypre 1D Matrix (size = m x 1)
    '''
    isComplex = MPI.COMM_WORLD.bcast(np.iscomplexobj(array), root=0)

    if isComplex:
       if array is None:
           rarray = None
           iarray = None
       else:
           rarray= array.real
           iarray= array.imag
       return  CHypreVec(Array2HypreVec(rarray, part),
                         Array2HypreVec(iarray, part))
    else:
       if array is None:
           rarray = None
       else:
           rarray= array
       return CHypreVec(Array2Hypre(rarray, part), None)

def CHypreVec2Array(array):
    from mpi4py import MPI
    myid     = MPI.COMM_WORLD.rank
    
    if array[0] is not None:
        r= HypreVec2Array(array[0])
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
            return r + 1j*i            
        else:
            return None
        
def CHypreMat2Coo(mat):
    print("CHYPREMat2Coo: deprecated,  Use class method !!!!")
    if mat.isComplex():
         return ToScipyCoo(mat.real) + 1j*ToScipyCoo(mat.imag)
    else:
         return ToScipyCoo(mat.real)

def LF2PyVec(rlf, ilf = None, horizontal= False):
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
        return CHypreVec(rv, iv, horizontal = horizontal)
    else:
        b1 = rlf.GetDataArray().copy()#; rlf.thisown = False
        if ilf is not None:
           b2 = ilf.GetDataArray()#; ilf.thisown = False
           b1 = b1+1j*b2
        if horizontal:
           return b1.reshape((1, -1))               
        else:
           return b1.reshape((-1, 1))
       
LinearForm2PyVector = LF2PyVec

def MfemVec2PyVec(rlf, ilf = None, horizontal= False):
    b1 = rlf.GetDataArray().copy()#; rlf.thisown = False
    if ilf is not None:
       b2 = ilf.GetDataArray()#; ilf.thisown = False
    else:
       b2 = None
       
    if MFEM_PAR:
        b1 = ToHypreParVec(b1)
        if b2 is not None:
           b2 = ToHypreParVec(b2)
        return CHypreVec(b1, b2, horizontal = horizontal)           
    else:
        if b2 is not None: b1 = b1+1j*b2
        if horizontal:
            return b1.reshape((1, -1))               
        else:
            return b1.reshape((-1, 1))
       
LinearForm2PyVector = LF2PyVec


def BF2PyMat(rbf, ibf = None, finalize = False):
    '''
    Convert pair of BilinearForms to CHypreMat or 
    ScipySparsematrix
    '''
    if finalize:
        rbf.Finalize()
        if ibf is not None: ibf.Finalize()
        
    if MFEM_PAR:
        M1 =  rbf.ParallelAssemble()
        M1.thisown = True        
        if ibf is not None:
            M2 =  ibf.ParallelAssemble()
            M2.thisown = True                    
        else:
            M2 = None
        return CHypreMat(M1,  M2)
    else:
        from mfem.common.sparse_utils import sparsemat_to_scipycsr        
        M1 =  rbf.SpMat()        
        if ibf is None:        
            return sparsemat_to_scipycsr(M1, dtype = float)
        if ibf is not None:
            M2 =  ibf.SpMat()
            m1 = sparsemat_to_scipycsr(M1, dtype = float).tolil()
            m2 = sparsemat_to_scipycsr(M2, dtype = complex).tolil()
            m = m1 + 1j*m2
            m = m.tocsr()
        return m

BilinearForm2PyMatix = BF2PyMat

def MfemMat2PyMat(M1, M2 = None):
    '''
    Convert pair of SpMat/HypreParCSR to CHypreMat or 
    ScipySparsematrix. This is simpler version of BF2PyMat, only 
    difference is it skippes convertion from BF to Matrix.
    '''
    from mfem.common.sparse_utils import sparsemat_to_scipycsr
    if MFEM_PAR:
        return CHypreMat(M1,  M2)
    else:
        if M2 is None:        
            return sparsemat_to_scipycsr(M1, dtype = float)
        else:
            m1 = sparsemat_to_scipycsr(M1, dtype = float).tolil()
            m2 = sparsemat_to_scipycsr(M2, dtype = complex).tolil()
            m = m1 + 1j*m2
            m = m.tocsr()
            return m




