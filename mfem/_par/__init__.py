import mfem

if mfem.mfem_mode is None:
   mfem.mfem_mode = 'parallel'
if mfem.mfem_mode == 'serial':
   raise ImportError("MFEM serial mode is already loaded")
debug_print = mfem.debug_print

from mpi4py import MPI

