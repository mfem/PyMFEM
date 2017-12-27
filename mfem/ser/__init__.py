import mfem

if mfem.mfem_mode is None:
   mfem.mfem_mode = 'serial'
if mfem.mfem_mode == 'parallel':
   raise ImportError("MFEM parallel mode is already loaded")

import sys, ctypes
## libmfem.a is linked only with _array.so
## this make sure that symbols are resovled
rtld_now = sys.getdlopenflags()
sys.setdlopenflags(ctypes.RTLD_GLOBAL|sys.getdlopenflags())

from  array import *
from  common_functions import *
from  socketstream import *
from  operators import *
from  blockoperator import *
from  blockvector import *
from  blockmatrix import *
from  coefficient import *
from  lininteg import *
from  mesh import *
from  fe_coll import *
from  vector import *
from  fespace import *
from  linearform import *
from  bilininteg import *
from  gridfunc import *
from  intrules import *
from  fe import *
from  ode import *
from  bilinearform import *
from  estimators import *
from  mesh_operators import *
from  sparsemat import *
from  densemat import *
from  solvers import *
from  sparsesmoothers import *
from  eltrans import *
from  geom import *
from  vertex import *
from  table import *
from  element import *
from  nonlininteg import *
from  nonlinearform import *

sys.setdlopenflags(rtld_now)

