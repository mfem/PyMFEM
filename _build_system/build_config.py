"""
Helper functions for setup.py
"""

import os
import sys
import configparser
from urllib import request
import itertools
import site
import re
import subprocess
import multiprocessing
import ssl
import tarfile
import shutil
from collections import namedtuple
from shutil import which as find_command

__all__ = ["print_config",
           "initialize_cmd_options",
           "cmd_options",
           "process_cmd_options",
           "configure_build",
           "clean_dist_info",
           ]

from build_utils import *
from build_consts import *
import build_globals as bglb


def print_config():
    print("----configuration----")
    print(" prefix", bglb.prefix)
    print(" when needed, the dependency (mfem/hypre/metis) will be installed under " +
          bglb.ext_prefix)
    print(" build mfem : " + ("Yes" if bglb.build_mfem else "No"))
    print(" build metis : " + ("Yes" if bglb.build_metis else "No"))
    print(" build hypre : " + ("Yes" if bglb.build_hypre else "No"))
    print(" build libceed : " + ("Yes" if bglb.build_libceed else "No"))
    print(" build gslib : " + ("Yes" if bglb.build_gslib else "No"))
    print(" call SWIG wrapper generator: " +
          ("Yes" if bglb.run_swig else "No"))
    print(" build serial wrapper: " + ("Yes" if bglb.build_serial else "No"))
    print(" build parallel wrapper : " +
          ("Yes" if bglb.build_parallel else "No"))

    print(" hypre prefix", bglb.hypre_prefix)
    print(" metis prefix", bglb.metis_prefix)
    print(" c compiler : " + bglb.cc_command)
    print(" c++ compiler : " + bglb.cxx_command)
    print(" mpi-c compiler : " + bglb.mpicc_command)
    print(" mpi-c++ compiler : " + bglb.mpicxx_command)

    print(" verbose : " + ("Yes" if bglb.verbose else "No"))
    print(" SWIG : " + swig_command)

    if bglb.blas_libraries != "":
        print(" BLAS libraries : " + bglb.blas_libraries)
    if bglb.lapack_libraries != "":
        print(" Lapack libraries : " + bglb.lapack_libraries)

    print("")


def clean_dist_info(wheeldir):
    if not os.path.isdir(wheeldir):
        return
    for x in os.listdir(wheeldir):
        if x.endswith(".dist-info"):
            fname = os.path.join(wheeldir, x)
            print("!!! removing existing ", fname)
            shutil.rmtree(fname)


def initialize_cmd_options(command_obj):
    command_obj.swig = False
    command_obj.skip_swig = False
    command_obj.ext_only = False

    command_obj.git_sshclone = False
    command_obj.skip_ext = False
    command_obj.with_parallel = False
    command_obj.no_serial = False
    command_obj.mfem_prefix = ''
    command_obj.mfems_prefix = ''
    command_obj.mfemp_prefix = ''
    command_obj.mfem_source = bglb.mfem_source
    command_obj.mfem_branch = ''
    command_obj.mfem_debug = False
    command_obj.mfem_build_miniapps = False
    command_obj.metis_prefix = ''
    command_obj.hypre_prefix = ''

    command_obj.with_cuda = False
    command_obj.with_cuda_hypre = False
    command_obj.cuda_arch = None
    command_obj.with_metis64 = False

    command_obj.with_pumi = False
    command_obj.pumi_prefix = ''

    command_obj.with_strumpack = False
    command_obj.strumpack_prefix = ''

    command_obj.with_suitesparse = False
    command_obj.suitesparse_prefix = ''

    command_obj.with_lapack = False
    command_obj.blas_libraries = ""
    command_obj.lapack_libraries = ""

    command_obj.with_libceed = False
    command_obj.libceed_prefix = ''
    command_obj.libceed_only = False

    command_obj.with_gslib = False
    command_obj.gslib_prefix = ''
    command_obj.gslib_only = False

    command_obj.CC = ''
    command_obj.CXX = ''
    command_obj.MPICC = ''
    command_obj.MPICXX = ''
    command_obj.vv = False

    command_obj.unverifiedSSL = False


cmd_options = [
    ('vv', None, 'More verbose output (CMAKE_VERBOSE_MAKEFILE etc)'),
    ('with-parallel', None, 'Installed both serial and parallel version'),
    ('no-serial', None, 'Skip building the serial wrapper'),
    ('mfem-prefix=', None, 'Specify locaiton of mfem' +
     'libmfem.so must exits under <mfem-prefix>/lib. ' +
     'This mode uses clean-swig + run-swig, unless mfem-prefix-no-swig is on'),
    ('mfemp-prefix=', None, 'Specify locaiton of parallel mfem ' +
     'libmfem.so must exits under <mfemp-prefix>/lib. ' +
     'Need to use it with mfem-prefix'),
    ('mfems-prefix=', None, 'Specify locaiton of serial mfem ' +
     'libmfem.so must exits under <mfems-prefix>/lib. ' +
     'Need to use it with mfem-prefix'),
    ('mfem-branch=', None, 'Specify branch of mfem' +
     'MFEM is cloned and built using the specfied branch '),
    ('mfem-source=', None, 'Specify mfem source location' +
     'MFEM source directory. Required to run-swig '),
    ('mfem-debug', None, 'Build MFME with MFEM_DEBUG enabled'),
    ('mfem-build-miniapps', None, 'build MFME Miniapps'),
    ('hypre-prefix=', None, 'Specify locaiton of hypre' +
     'libHYPRE.so must exits under <hypre-prefix>/lib'),
    ('metis-prefix=', None, 'Specify locaiton of metis' +
     'libmetis.so must exits under <metis-prefix>/lib'),
    ('git-sshclone', None, 'Use SSH for git clone' +
     'try if default git clone using https fails (need Github account and setting for SSH)'),
    ('swig', None, 'Run Swig and exit'),
    ('skip-swig', None,
     'Skip running swig (used when wrapper is generated for the MFEM C++ library to be used'),
    ('ext-only', None, 'Build metis, hypre, mfem(C++) only'),
    ('skip-ext', None, 'Skip building metis, hypre, mfem(C++) only'),
    ('CC=', None, 'c compiler'),
    ('CXX=', None, 'c++ compiler'),
    ('MPICC=', None, 'mpic compiler'),
    ('MPICXX=', None, 'mpic++ compiler'),
    ('unverifiedSSL', None, 'use unverified SSL context for downloading'),
    ('with-cuda', None, 'enable cuda'),
    ('with-cuda-hypre', None, 'enable cuda in hypre'),
    ('cuda-arch=', None, 'set cuda compute capability. Ex if A100, set to 80'),
    ('with-metis64', None, 'use 64bit int in metis'),
    ('with-pumi', None, 'enable pumi (parallel only)'),
    ('pumi-prefix=', None, 'Specify locaiton of pumi'),
    ('with-suitesparse', None,
     'build MFEM with suitesparse (MFEM_USE_SUITESPARSE=YES) (parallel only)'),
    ('suitesparse-prefix=', None,
     'Specify locaiton of suitesparse (=SuiteSparse_DIR)'),
    ('with-libceed', None, 'enable libceed'),
    ('libceed-prefix=', None, 'Specify locaiton of libceed'),
    ('libceed-only', None, 'Build libceed only'),
    ('gslib-prefix=', None, 'Specify locaiton of gslib'),
    ('with-gslib', None, 'enable gslib'),
    ('gslib-only', None, 'Build gslib only'),
    ('with-strumpack', None, 'enable strumpack (parallel only)'),
    ('strumpack-prefix=', None, 'Specify locaiton of strumpack'),
    ('with-lapack', None, 'build MFEM with lapack'),
    ('blas-libraries=', None, 'Specify locaiton of Blas library (used to build MFEM)'),
    ('lapack-libraries=', None,
     'Specify locaiton of Lapack library (used to build MFEM)'),
]


def process_cmd_options(command_obj, cfs):
    '''
    called when install workflow is used
    '''
    cc = cfs.pop("CC", "")
    if cc != "":
        command_obj.cc_command = cc

    cc = cfs.pop("CXX", "")
    if cc != "":
        command_obj.cxx_command = cc

    cc = cfs.pop("MPICC", "")
    if cc != "":
        command_obj.mpicc_command = cc

    cc = cfs.pop("MPICXX", "")
    if cc != "":
        command_obj.mpicxx_command = cc

    for item in cmd_options:
        param, _none, hit = item
        attr = "_".join(param.split("-"))

        if param.endswith("="):
            param = param[:-1]
            attr = attr[:-1]
            value = cfs.pop(param, "")
            if value != "":
                if not hasattr(command_obj, attr):
                    assert False, str(command_obj) + " does not have " + attr
                setattr(command_obj, attr, value)
        else:
            value = cfs.pop(param, "No")
            if not hasattr(command_obj, attr):
                assert False, str(command_obj) + " does not have " + attr

            if value.upper() in ("YES", "TRUE", "1"):
                setattr(command_obj, attr, True)
            else:
                setattr(command_obj, attr, False)


def process_setup_options(command_obj, args):
    for item in args:
        if item.startswith('--'):
            item = item[2:]
        if item.startswith('-'):
            item = item[1:]

        if len(item.split('=')) == 2:
            param = item.split('=')[0]
            value = item.split('=')[1]
        else:
            param = item.strip()
            value = True
        attr = "_".join(param.split("-"))

        setattr(command_obj, attr, value)


def configure_install(self):
    '''
    called when install workflow is used

    '''
    if sys.argv[0] == 'setup.py' and sys.argv[1] == 'install':
        process_setup_options(self, sys.argv[2:])
    else:
        if bglb.verbose:
            print("!!!!!!!!  command-line input (pip): ", bglb.cfs)
        process_cmd_options(self, bglb.cfs)

    bglb.verbose = bool(self.vv) if not bglb.verbose else bglb.verbose
    if bglb.dry_run:
        bglb.verbose = True

    bglb.git_sshclone = bool(self.git_sshclone)

    bglb.mfem_source = abspath(self.mfem_source)

    bglb.skip_ext = bool(self.skip_ext)
    bglb.skip_swig = bool(self.skip_swig)

    bglb.swig_only = bool(self.swig)
    bglb.ext_only = bool(self.ext_only)

    bglb.metis_64 = bool(self.with_metis64)
    bglb.enable_pumi = bool(self.with_pumi)
    bglb.enable_strumpack = bool(self.with_strumpack)
    bglb.enable_cuda = bool(self.with_cuda)
    bglb.enable_cuda_hypre = bool(self.with_cuda_hypre)
    if self.cuda_arch is not None:
        bglb.cuda_arch = self.cuda_arch
    bglb.enable_libceed = bool(self.with_libceed)
    bglb.libceed_only = bool(self.libceed_only)
    bglb.enable_gslib = bool(self.with_gslib)
    bglb.gslib_only = bool(self.gslib_only)
    bglb.enable_suitesparse = bool(self.with_suitesparse)
    bglb.enable_lapack = bool(self.with_lapack)

    # controlls PyMFEM parallel
    bglb.build_parallel = bool(self.with_parallel)
    bglb.build_serial = not bool(self.no_serial)

    bglb.clean_swig = True
    bglb.run_swig = True
    bglb.run_swig_parallel = bool(self.with_parallel)

    bglb.mfem_debug = bool(self.mfem_debug)
    bglb.mfem_build_miniapps = bool(self.mfem_build_miniapps)

    if bglb.build_serial:
        bglb.build_serial = (not bglb.swig_only and not bglb.ext_only)

    if bglb.build_parallel:
        try:
            import mpi4py
        except ImportError:
            assert False, "Can not import mpi4py"

    if self.mfem_prefix != '':
        bglb.mfem_prefix = abspath(self.mfem_prefix)
        bglb.mfems_prefix = abspath(self.mfem_prefix)
        bglb.mfemp_prefix = abspath(self.mfem_prefix)
        if self.mfems_prefix != '':
            bglb.mfems_prefix = abspath(self.mfems_prefix)
        if self.mfemp_prefix != '':
            bglb.mfemp_prefix = abspath(self.mfemp_prefix)

        check = find_libpath_from_prefix('mfem', bglb.mfems_prefix)
        assert check != '', "libmfem.so is not found in the specified <path>/lib"
        check = find_libpath_from_prefix('mfem', bglb.mfemp_prefix)
        assert check != '', "libmfem.so is not found in the specified <path>/lib"

        bglb.mfem_outside = True
        bglb.build_mfem = False
        hypre_prefix = bglb.mfem_prefix
        metis_prefix = bglb.mfem_prefix

        if bglb.swig_only:
            bglb.clean_swig = False

    else:
        bglb.mfem_outside = False
        bglb.build_mfem = True
        bglb.build_mfemp = bglb.build_parallel
        bglb.build_hypre = bglb.build_parallel
        bglb.build_metis = bglb.build_parallel or bglb.enable_suitesparse

        if bglb.ext_prefix == '':
            bglb.ext_prefix = external_install_prefix(bglb.prefix)
        bglb.hypre_prefix = os.path.join(bglb.ext_prefix)
        bglb.metis_prefix = os.path.join(bglb.ext_prefix)

        bglb.mfem_prefix = bglb.ext_prefix
        bglb.mfems_prefix = os.path.join(bglb.ext_prefix, 'ser')
        bglb.mfemp_prefix = os.path.join(bglb.ext_prefix, 'par')
        # enable_gslib = True

    if self.mfem_branch != '':
        bglb.mfem_branch = self.mfem_branch

    if self.hypre_prefix != '':
        check = find_libpath_from_prefix('HYPRE', self.hypre_prefix)
        assert check != '', "libHYPRE.so is not found in the specified <path>/lib or lib64"
        hypre_prefix = os.path.expanduser(self.hypre_prefix)
        build_hypre = False

    if self.metis_prefix != '':
        check = find_libpath_from_prefix('metis', self.metis_prefix)
        assert check != '', "libmetis.so is not found in the specified <path>/lib or lib64"
        bglb.metis_prefix = os.path.expanduser(self.metis_prefix)
        bglb.build_metis = False

    if bglb.enable_libceed or bglb.libceed_only:
        if self.libceed_prefix != '':
            bglb.libceed_prefix = os.path.expanduser(self.libceed_prefix)
            bglb.build_libceed = False
        else:
            bglb.libceed_prefix = bglb.mfem_prefix
            bglb.build_libceed = True

    if bglb.enable_gslib or bglb.gslib_only:
        if self.gslib_prefix != '':
            bglb.build_gslib = False
            bglb.gslibs_prefix = os.path.expanduser(self.gslib_prefix)
            bglb.gslibp_prefix = os.path.expanduser(self.gslib_prefix)
        else:
            bglb.gslibs_prefix = bglb.mfems_prefix
            bglb.gslibp_prefix = bglb.mfemp_prefix
            bglb.build_gslib = True

    if bglb.enable_suitesparse and self.suitesparse_prefix != '':
        bglb.suitesparse_prefix = self.suitesparse_prefix

    if self.pumi_prefix != '':
        bglb.pumi_prefix = abspath(self.pumi_prefix)
    else:
        bglb.pumi_prefix = bglb.mfem_prefix

    if self.strumpack_prefix != '':
        bglb.strumpack_prefix = abspath(self.strumpack_prefix)
    else:
        bglb.strumpack_prefix = bglb.mfem_prefix

    if bglb.enable_cuda:
        nvcc = find_command('nvcc')
        bglb.cuda_prefix = os.path.dirname(os.path.dirname(nvcc))

    if self.CC != '':
        bglb.cc_command = self.CC
    if self.CXX != '':
        bglb.cxx_command = self.CXX
    if self.MPICC != '':
        bglb.mpicc_command = self.MPICC
    if self.MPICXX != '':
        bglb.mpicxx_command = self.MPICXX

    if self.blas_libraries != "":
        bglb.blas_libraries = self.blas_libraries
    if self.lapack_libraries != "":
        bglb.lapack_libraries = self.lapack_libraries

    if bglb.swig_only:
        bglb.build_serial = False
        bglb.build_parallel = False
        bglb.clean_swig = False
        bglb.keep_temp = True
        bglb.skip_ext = True

    if bglb.skip_ext:
        bglb.build_metis = False
        bglb.build_hypre = False
        bglb.build_mfem = False
        bglb.build_mfemp = False
        bglb.build_libceed = False
        bglb.build_gslib = False

    if bglb.skip_swig:
        bglb.clean_swig = False
        bglb.run_swig = False

    if bglb.ext_only:
        bglb.clean_swig = False
        bglb.run_swig = False
        bglb.build_serial = False
        bglb.build_parallel = False
        bglb.keep_temp = True


    if bglb.libceed_only:
        bglb.clean_swig = False
        bglb.run_swig = False
        bglb.build_mfem = False
        bglb.build_mfemp = False
        bglb.build_metis = False
        bglb.build_hypre = False
        bglb.build_gslib = False
        bglb.build_serial = False
        bglb.build_parallel = False
        bglb.build_libceed = True
        bglb.keep_temp = True


    if bglb.gslib_only:
        bglb.clean_swig = False
        bglb.run_swig = False
        bglb.build_mfem = False
        bglb.build_mfemp = False
        bglb.build_metis = False
        bglb.build_hypre = False
        bglb.build_serial = False
        bglb.build_libceed = False
        bglb.build_gslib = True
        bglb.keep_temp = True


    bglb.is_configured = True

configure_build = configure_install
