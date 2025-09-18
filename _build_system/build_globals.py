# ----------------------------------------------------------------------------------------
# Global build parameters
# ----------------------------------------------------------------------------------------
import os

is_configured = False
prefix = ''

verbose = True
git_sshclone = False
swig_only = False
run_swig = False
clean_swig = False
build_mfem = False
mfem_branch = None
build_mfemp = False
build_metis = False
build_hypre = False
build_libceed = False
build_gslib = False
build_parallel = False
build_serial = False

ext_prefix = ''
mfem_outside = False
mfems_prefix = ''
mfemp_prefix = ''
mfem_source = os.path.join(os.path.dirname(__file__), "..", "external", "mfem")

metis_prefix = ''
hypre_prefix = ''

enable_cuda = False
enable_cuda_hypre = False
cuda_prefix = ''
cuda_arch = ''
enable_pumi = False
pumi_prefix = ''
enable_strumpack = False
strumpack_prefix = ''
enable_libceed = False
libceed_prefix = ''
libceed_only = False
enable_gslib = False
gslibs_prefix = ''
gslibp_prefix = ''
gslib_only = False
mfem_debug = False
mfem_build_miniapps = True

enable_suitesparse = False
suitesparse_prefix = "/usr/"

enable_lapack = False
blas_libraries = ""
lapack_libraries = ""

dry_run = False
do_bdist_wheel = False
bdist_wheel_dir = ''

keep_temp = False

use_unverifed_SSL = False if os.getenv(
    "unverifedSSL") is None else os.getenv("unverifiedSSL")

use_metis_gklib = False
metis_64 = False

# ----------------------------------------------------------------------------------------
#   command line configuration parameter (pip -C)
# ----------------------------------------------------------------------------------------

cfs = {}

# ----------------------------------------------------------------------------------------
#   enviromental variables.
# ----------------------------------------------------------------------------------------
cc_command = 'cc' if os.getenv("CC") is None else os.getenv("CC")
cxx_command = 'c++' if os.getenv("CC") is None else os.getenv("CXX")
mpicc_command = 'mpicc' if os.getenv("MPICC") is None else os.getenv("MPICC")
mpicxx_command = 'mpic++' if os.getenv(
    "MPICXX") is None else os.getenv("MPICXX")
cxxstd_flag = '-std=c++17' if os.getenv(
    "CXXSTDFLAG") is None else os.getenv("CXXSTDFLAG")

# location of MPI.h. Usually it is not needed, as long as MPI compiler can be used
#
mpiinc = '' if os.getenv("MPIINC") is None else os.getenv("MPIINC")
