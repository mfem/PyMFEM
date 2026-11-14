import mfem

if mfem.mfem_mode == 'serial':
   raise ImportError("MFEM serial mode is already loaded")
debug_print = mfem.debug_print

from mfem.common.check_parallel import check_parallel_version
from mfem.common.find_libmpi import check_mpi_runtime

check_parallel_version()
check_mpi_runtime()

from mpi4py import MPI
mfem.mfem_mode = 'parallel'
