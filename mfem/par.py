import mfem

if mfem.mfem_mode is None:
   mfem.mfem_mode = 'parallel'
if mfem.mfem_mode == 'serial':
   raise ImportError("MFEM serial mode is already loaded")
debug_print = mfem.debug_print

from mpi4py import MPI

from  mfem._par.cpointers import *
from  mfem._par.globals import *
from  mfem._par.mem_manager import *
from  mfem._par.device import *
from  mfem._par.hash import *
from  mfem._par.point import *
from  mfem._par.segment import *
from  mfem._par.triangle import *
from  mfem._par.quadrilateral import *
from  mfem._par.wedge import *
from  mfem._par.tetrahedron import *
from  mfem._par.hexahedron import *
from  mfem._par.array import *
from  mfem._par.common_functions import *
from  mfem._par.socketstream import *
from  mfem._par.operators import *
from  mfem._par.blockoperator import *
from  mfem._par.blockvector import *
from  mfem._par.blockmatrix import *
from  mfem._par.coefficient import *
from  mfem._par.lininteg import *
from  mfem._par.handle import *
from  mfem._par.mesh import *
from  mfem._par.fe_coll import *
from  mfem._par.vector import *
from  mfem._par.complex_operator import *
from  mfem._par.complex_fem import *
from  mfem._par.fespace import *
from  mfem._par.linearform import *
from  mfem._par.bilininteg import *
from  mfem._par.gridfunc import *
from  mfem._par.intrules import *
from  mfem._par.fe import *
from  mfem._par.ode import *
from  mfem._par.bilinearform import *
from  mfem._par.estimators import *
from  mfem._par.mesh_operators import *
from  mfem._par.sparsemat import *
from  mfem._par.densemat import *
from  mfem._par.solvers import *
from  mfem._par.sparsesmoothers import *
from  mfem._par.eltrans import *
from  mfem._par.geom import *
from  mfem._par.vertex import *
from  mfem._par.table import *
from  mfem._par.element import *
from  mfem._par.nonlininteg import *
from  mfem._par.nonlinearform import *
from  mfem._par.ncmesh import *
from  mfem._par.pmesh import *
from  mfem._par.pfespace import *
from  mfem._par.plinearform import *
from  mfem._par.pbilinearform import *
from  mfem._par.pnonlinearform import *
from  mfem._par.pgridfunc import *
from  mfem._par.hypre import *
from  mfem._par.stable3d import *
from  mfem._par.vtk import *
from  mfem._par.datacollection import *
from  mfem._par.io_stream import wFILE, STDOUT
from  mfem._par.fespacehierarchy import *
from  mfem._par.multigrid import *
from  mfem._par.constraints import *
from  mfem._par.transfer import *

from  mfem._par.fe_base import *
from  mfem._par.fe_h1 import *
from  mfem._par.fe_l2 import *
from  mfem._par.fe_nd import *
from  mfem._par.fe_rt import *
from  mfem._par.fe_ser import *
from  mfem._par.fe_fixed_order import *
from  mfem._par.fe_pos import *
from  mfem._par.fe_nurbs import *
from  mfem._par.doftrans import *
from  mfem._par.std_vectors import *

try:
   from  mfem._par.gslib import *
except:
   pass

import mfem._par.array as array
import mfem._par.blockoperator as blockoperator
import mfem._par.coefficient as coefficient
import mfem._par.cpointers as cpointers
import mfem._par.densemat as densemat
import mfem._par.error as error
import mfem._par.fe as fe
import mfem._par.fe_coll as fe_coll
import mfem._par.fespace as fespace
import mfem._par.geom as geom
import mfem._par.gridfunc as gridfunc
import mfem._par.hypre as hypre
import mfem._par.intrules as intrules
import mfem._par.mesh as mesh
import mfem._par.pgridfunc as pgridfunc
import mfem._par.solvers as solvers
import mfem._par.vector as vector
import mfem._par.sparsemat as sparsemat

import mfem._par.tmop_modules as tmop



try:
    import mfem._par.dist_solver as dist_solver
except:
    pass

try:   
    from mfem._par.schwarz import (SchwarzSmoother,
                                   ComplexSchwarzSmoother)
except:
    pass

try:
   import mfem._par.pumi as pumi
   from mfem._par.pumi import *
except ImportError:
   pass
