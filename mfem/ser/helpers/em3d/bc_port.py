'''
define TE mode shape
'''
import mfem
class TE_Mode(mfem.VectorPyCoefficient):
   def __init__(self, c0, norm,  v1, v2, m = 0, n = 1, field='E'):
       self.m = m
       self.n = 0
       self.c0 = c0      # center (x, y, z) of port
       self.norm = norm  # normal direction of foward wave
       self.v1 = v1      # direction of v1 (longer edge)
       self.v2 = v2      # direction of v2 (shorter edge)
       self.field = 'E'
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi), 0, 0]
   
class TE_H(mfem.VectorPyCoefficient):
   def __init__(self, m = 0, n = 1):
       self.m = m
       self.n = 0
   def EvalValue(self, x):
       return [0., 0., k/mu0*np.sin(x[2]/b*np.pi)]       

class Port():
    def __init__(self, kbdr=-1, mode = 'te01', kport=1, ext=1):
        '''
        mode: te_m_n
              tm_m_n
              tem_m_n
              c_m_n

        ext: !=0: excitation
               0: passive
        '''
        self.kbdr = kbdr    
        self.mode = mode
        self.kport = kport
        self.ext = ext

    def CollectMeshInfo(self, mesh):

    def Assemble(fes):
        '''
        assemble bilinear and linear form
        [A    c1] [x]    [b]
        [       ] [ ]  = [ ]
        [c2   c3] [l]    [m]
        '''
        '''

