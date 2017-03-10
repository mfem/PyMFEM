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
from scipy.sparse  import csr_matrix
from mfem.common.parcsr_extra import *


from mpi4py import MPI
myid     = MPI.COMM_WORLD.rank
class CHypreVec(list):
    def __init__(self, r = None, i = None):
        list.__init__(self, [None]*2)
        if isinstance(r, np.ndarray):
            self[0] = ToHypreParVec(r)
        else:
            self[0] = r
        if isinstance(i, np.ndarray):
            self[1] = ToHypreParVec(i)
        else:
            self[1] = i
            
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
        return 1, self[0].GlobalSize()
            
    def isComplex(self):
        return not (self[1] is None)
            
    def __imul__(self, other):
        if self[0] is not None:
            self[0] *= other
        if self[1] is not None:
            self[1] *= other
            
    def dot(self, other):
        if not isinstance(other, CHypreVec):
             raise ValueError(
                   "argument should be CHypreVec")
        return InnerProductComplex(A, B)
    
    def set_element(self, i, v):
        if self[0] is not None:
            part = self[0].GetPartitioningArray()
        else:
            part = self[1].GetPartitioningArray()

        if part[0] <= i and i < part[1]:
            v = complex(v)
            if self[0] is not None:
                self[0][int(i - part[0])] = v.real
            if self[1] is not None:            
                self[1][int(i - part[0])] = v.imag
                
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

    def isComplex(self):
        return not (self[1] is None)

    def __mul__(self, other): # A * B or A * v
        if isinstance(other, CHypreMat):
            return CHypreMat(*ParMultComplex(self, other))
        if isinstance(other, CHypreVec):
            return CHypreVec(*ParMultVecComplex(self, other))
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
         
        if self[0] is not None and other[0] is None:
            r = ToHypreParCSR((ToScipyCoo(self[0])+ ToScipyCoo(other[0])).tocsr())            
        elif self[0] is not None:
            r = self[0]
        elif other[0] is not None:
            r = other[0]            
        else:
            r = None

        if self[1] is not None and other[1] is None:
            i = ToHypreParCSR((ToScipyCoo(self[1])+ ToScipyCoo(other[1])).tocsr())            
        elif self[1] is not None:
            i = self[1]
        elif other[1] is not None:
            i = other[1]            
        else:
            i = None
            
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
        #ret = self*B
        #return B.transpose().conj() * ret
        ret = B.conj().transpose()*self
        ret = ret * B.conj()  # this should put B back to original
        return ret

    def setDiag(self, idx, value):
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

    @property
    def nnz(self):
        if self[0] is not None and self[1] is not None:
            return self[0].NNZ(), self[1].NNZ()
        if self[0] is not None: return self[0].NNZ()
        if self[1] is not None: return self[1].NNZ()
    
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
        return self.m(), self.n()


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
    if mat.isComplex():
         return ToScipyCoo(mat.real) + 1j*ToScipyCoo(mat.imag)
    else:
         return ToScipyCoo(mat.real)

       
