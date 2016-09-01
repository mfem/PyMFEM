'''
  CHypre (Complex Hypre)
  
  container object to support complex using
  real value hypre

  it should work with pure real or pure imaginary
  case too.
'''
import numpy as np
from scipy.sparse  import csr_matrix
from mfem.common.parcsr_extra import *


#ParMultComplex, RapComplex, Array2Hypre, Hypre2Array,
#TransposeComplex, ResetHypreDiag, ToScipyCoo,
#get_row_partitioning, ParAdd, get_col_partitioning, Conj

from mpi4py import MPI
myid     = MPI.COMM_WORLD.rank

class CHypreMat(list):
    def __init__(self, r = None, i = None):
        list.__init__(self, [None]*2)
        if isinstance(r, csr_matrix):
            self[0] = ToHypreParCSR(r)
        else:
            self[0] = r
        if isinstance(i, csr_matrix):
            self[1] = ToHypreParCSR(i)
        else:
            self[1] = i

    def isComplex(self):
        return not (self[1] is None)

    def __mul__(self, other): # A * B
        if not isinstance(other, CHypreMat):
             raise ValueError(
                   "argument should be CHypreMat")
        return CHypreMat(*ParMultComplex(self, other))
    
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

    def dot(self, other):
        return self.__mul__(other)
    
    def transpose(self):
        return CHypreMat(self[0].Transpose(), self[1].Transpose())

    def conj(self, inplace = True):
        '''
        complex conjugate
          if copy is on, imaginary part becomes different object
        ''' 
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


    def nnz(self):
        ans = []
        if self[0] is not None: ans.append(self[0].NNZ())
        if self[1] is not None: ans.append(self[1].NNZ())
        return ans

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
        
        #dprint3('NNZ ', M1.NNZ(), ' ' ,  M2.NNZ(), ' ',  M1.M(), ' ', M1.N())            

def Array2CHypre(array, part):
    isComplex = MPI.COMM_WORLD.bcast(np.iscomplexobj(array), root=0)

    if isComplex:
       if array is None:
           rarray = None
           iarray = None
       else:
           rarray= array.real
           iarray= array.imag
       return  CHypreMat(Array2Hypre(rarray, part),
                     Array2Hypre(iarray, part))
    else:
       if array is None:
           rarray = None
       else:
           rarray= array
       return CHypreMat(Array2Hypre(rarray, part), None)

def CHypre2Array(array):
    if array[0] is not None:
        r= Hypre2Array(array[0])
    else:
        r = 0.0
    if array[1] is None:
        return r
    else:
        i = Hypre2Array(array[1])
    if i is None:
        return r
    else:
        return r + 1j*i


       
