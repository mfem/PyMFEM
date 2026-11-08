"""
Global build parameters with fixed member names and mutable values.
"""
import os
from dataclasses import dataclass, field


@dataclass(slots=True)
class BuildGlobal:
    """
    Declare every build parameter before configuration resolves its value.
    """

    is_configured: bool = False
    prefix: str = ''

    verbose: bool = True
    git_sshclone: bool = False
    swig_only: bool = False
    run_swig: bool = False
    clean_swig: bool = False
    build_mfem: bool = False
    mfem_branch: str | None = None
    build_mfemp: bool = False
    build_metis: bool = False
    build_hypre: bool = False
    build_libceed: bool = False
    build_gslib: bool = False
    build_parallel: bool = False
    build_serial: bool = False

    ext_prefix: str = ''
    mfem_outside: bool = False
    mfems_prefix: str = ''
    mfemp_prefix: str = ''
    mfem_source: str = os.path.join(os.path.dirname(__file__), "..", "external", "mfem")

    metis_prefix: str = ''
    hypre_prefix: str = ''

    enable_cuda: bool = False
    enable_cuda_hypre: bool = False
    cuda_prefix: str = ''
    cuda_arch: str = ''
    enable_pumi: bool = False
    pumi_prefix: str = ''
    enable_strumpack: bool = False
    strumpack_prefix: str = ''
    enable_libceed: bool = False
    libceed_prefix: str = ''
    libceed_only: bool = False
    enable_gslib: bool = False
    gslibs_prefix: str = ''
    gslibp_prefix: str = ''
    gslib_only: bool = False
    mfem_debug: bool = False
    mfem_build_miniapps: bool = True

    enable_suitesparse: bool = False
    suitesparse_prefix: str = "/usr/"

    enable_lapack: bool = False
    blas_libraries: str = ""
    lapack_libraries: str = ""

    dry_run: bool = False
    do_bdist_wheel: bool = False
    bdist_wheel_dir: str = ''

    keep_temp: bool = False

    use_unverifed_SSL: bool | str | None = False if os.getenv(
        "unverifedSSL") is None else os.getenv("unverifiedSSL")

    use_metis_gklib: bool = False
    metis_64: bool = False

    # ----------------------------------------------------------------------------------------
    #   command line configuration parameter (pip -C)
    # ----------------------------------------------------------------------------------------

    cfs: dict = field(default_factory=dict)

    # ----------------------------------------------------------------------------------------
    #   enviromental variables.
    # ----------------------------------------------------------------------------------------
    cc_command: str | None = 'cc' if os.getenv("CC") is None else os.getenv("CC")
    cxx_command: str | None = 'c++' if os.getenv("CC") is None else os.getenv("CXX")
    mpicc_command: str | None = 'mpicc' if os.getenv("MPICC") is None else os.getenv("MPICC")
    mpicxx_command: str | None = 'mpic++' if os.getenv(
        "MPICXX") is None else os.getenv("MPICXX")
    cxxstd_flag: str = '-std=c++17' if os.getenv(
        "CXXSTDFLAG") is None else os.getenv("CXXSTDFLAG")

    # location of MPI.h. Usually it is not needed, as long as MPI compiler can be used
    #
    mpiinc: str = '' if os.getenv("MPIINC") is None else os.getenv("MPIINC")

    # Fields assigned by setup.py and configure_build().
    build_py_done: bool = False
    skip_ext: bool = False
    skip_swig: bool = False
    ext_only: bool = False
    run_swig_parallel: bool = False
    mfem_prefix: str = ''
    mfem_miniapps: bool = True


bglb = BuildGlobal()
