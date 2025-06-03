# ----------------------------------------------------------------------------------------
# Global build parameters
# ----------------------------------------------------------------------------------------
import os

is_configured = False
prefix = ''

verbose = True
git_sshclone = False
swig_only = False
skip_install = False
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

use_unverifed_SSL = False if os.getenv(
    "unverifedSSL") is None else os.getenv("unverifiedSSL")

use_metis_gklib = False
