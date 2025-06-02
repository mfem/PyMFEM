"""
  MFEM + PyMFEM (finite element method library)
"""

import sys
import os
import site
import re
import shutil
import subprocess
from shutil import which as find_command

from sys import platform
import multiprocessing
from multiprocessing import Pool
# force fork instead of spawn on MacOS to avoid race condition on mfem/__init__.py
if platform == "darwin":
    multiprocessing.set_start_method("fork")

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.install import install as _install
from setuptools.command.install_egg_info import install_egg_info as _install_egg_info
from setuptools.command.install_lib import install_lib as _install_lib
from setuptools.command.install_scripts import install_scripts as _install_scripts

from distutils.command.clean import clean as _clean

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_build_system"))

from build_utils import *
from build_consts import *
from build_pymfem import *
from build_hypre import *
from build_metis import *
from build_mfem import *

import build_globals as bglb

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------

from sys import platform
haveWheel = False

if platform == "linux" or platform == "linux2":
    try:
        from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
        haveWheel = True
    except ImportError:
        print("Skipping wheel build; wheel not installed.")

if haveWheel:

    class BdistWheel(_bdist_wheel):
        '''
        Wheel build performs SWIG + Serial in Default.
        --skip-build option skip building entirely.
        '''
        def finalize_options(self):
            def _has_ext_modules():
                return True
            self.distribution.has_ext_modules = _has_ext_modules
            _bdist_wheel.finalize_options(self)

        def run(self):
            print("!!!!! Entering BdistWheel::Run")
            import build_globals as bglb
            
            if not bglb.is_configured:
                print('running config')
                configure_bdist(self)
                print_config()
            self.run_command("build")
            _bdist_wheel.run(self)




def clean_so(all=None):
    python = sys.executable
    command = [python, "setup.py", "clean"]
    if all == 1:
        command.append("--all")

    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))
    for f in os.listdir():
        if f.endswith('.so'):
            os.remove(f)
        if f.endswith('.dylib'):
            os.remove(f)
    make_call(command)

    chdir(os.path.join(rootdir, 'mfem', '_par'))
    for f in os.listdir():
        if f.endswith('.so'):
            os.remove(f)
        if f.endswith('.dylib'):
            os.remove(f)
    make_call(command)

    chdir(pwd)


def print_config():
    import build_globals as bglb    

    print("----configuration----")
    print(" prefix", bglb.prefix)
    print(" when needed, the dependency (mfem/hypre/metis) will be installed under " +
          bglb.ext_prefix)
    print(" build mfem : " + ("Yes" if bglb.build_mfem else "No"))
    print(" build metis : " + ("Yes" if bglb.build_metis else "No"))
    print(" build hypre : " + ("Yes" if bglb.build_hypre else "No"))
    print(" build libceed : " + ("Yes" if bglb.build_libceed else "No"))
    print(" build gslib : " + ("Yes" if bglb.build_gslib else "No"))
    print(" call SWIG wrapper generator: " + ("Yes" if bglb.run_swig else "No"))
    print(" build serial wrapper: " + ("Yes" if bglb.build_serial else "No"))
    print(" build parallel wrapper : " + ("Yes" if bglb.build_parallel else "No"))

    print(" hypre prefix", bglb.hypre_prefix)
    print(" metis prefix", bglb.metis_prefix)
    print(" c compiler : " + cc_command)
    print(" c++ compiler : " + cxx_command)
    print(" mpi-c compiler : " + mpicc_command)
    print(" mpi-c++ compiler : " + mpicxx_command)

    print(" verbose : " + ("Yes" if bglb.verbose else "No"))
    print(" SWIG : " + swig_command)

    if bglb.blas_libraries != "":
        print(" BLAS libraries : " + bglb.blas_libraries)
    if bglb.lapack_libraries != "":
        print(" Lapack libraries : " + bglb.lapack_libraries)

    print("")


def configure_install(self):
    '''
    called when install workflow is used
    '''
    verbose = bool(self.vv) if bglb.verbose == -1 else bglb.verbose
    dry_run = bool(self.dry_run) if bglb.dry_run == -1 else bglb.dry_run
    if bglb.dry_run:
        bglb.verbose = True

    bglb.git_sshclone = bool(self.git_sshclone)

    bglb.prefix = abspath(self.prefix)
    bglb.mfem_source = abspath(self.mfem_source)

    bglb.skip_ext = bool(self.skip_ext)
    bglb.skip_install = bool(self.build_only)
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

    bglb.build_parallel = bool(self.with_parallel)     # controlls PyMFEM parallel
    bglb.build_serial = not bool(self.no_serial)

    bglb.clean_swig = True
    bglb.run_swig = True

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

        bglb.build_mfem = False
        hypre_prefix = bglb.mfem_prefix
        metis_prefix = bglb.mfem_prefix

        if bglb.swig_only:
            bglb.clean_swig = False

    else:
        bglb.build_mfem = True
        bglb.build_mfemp = bglb.build_parallel
        bglb.build_hypre = bglb.build_parallel
        bglb.build_metis = bglb.build_parallel or bglb.enable_suitesparse

        print("ext_prefix", bglb.ext_prefix)
        if bglb.ext_prefix == '':
            bglb.ext_prefix = external_install_prefix(prefix)
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
        bglb.pumi_prefix = mfem_prefix

    if self.strumpack_prefix != '':
        bglb.strumpack_prefix = abspath(self.strumpack_prefix)
    else:
        bglb.strumpack_prefix = bglb.mfem_prefix

    if bglb.enable_cuda:
        nvcc = find_command('nvcc')
        bglb.cuda_prefix = os.path.dirname(os.path.dirname(nvcc))

    '''
    this has to be handled differently
    if self.CC != '':
        cc_command = self.CC
    if self.CXX != '':
        cxx_command = self.CXX
    if self.MPICC != '':
        mpicc_command = self.MPICC
    if self.MPICXX != '':
        mpicxx_command = self.MPICXX
    '''

    if self.blas_libraries != "":
        bglb.blas_libraries = self.blas_libraries
    if self.lapack_libraries != "":
        bglb.lapack_libraries = self.lapack_libraries

    if skip_ext:
        bglb.build_metis = False
        bglb.build_hypre = False
        bglb.build_mfem = False
        bglb.build_mfemp = False
        bglb.build_libceed = False
        bglb.build_gslib = False

    if self.skip_swig:
        bglb.clean_swig = False
        bglb.run_swig = False

    if bglb.swig_only:
        bglb.build_serial = False
        bglb.clean_swig = False

    if bglb.ext_only:
        bglb.clean_swig = False
        bglb.run_swig = False
        bglb.build_serial = False
        bglb.build_parallel = False
        bglb.skip_install = True

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
        bglb.skip_install = True

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
        bglb.skip_install = True

    bglb.is_configured = True


def configure_bdist(self):
    '''
    called when bdist workflow is used
    '''
    
    bglb.dry_run = bool(self.dry_run) if bglb.dry_run == -1 else bglb.dry_run
    bglb.verbose = bool(self.verbose) if bglb.verbose == -1 else bglb.verbose

    bglb.prefix = abspath(self.bdist_dir)
    bglb.prefix = abspath(self.bdist_dir)

    bglb.build_parallel = False

    if self.skip_build == 1:
        bglb.build_mfem = False
        bglb.build_serial = False
        bglb.run_swig = False
    else:
        bglb.build_mfem = True
        bglb.build_serial = True
        # build_gslib = True
        bglb.run_swig = True

    bglb.is_configured = True
    bglb.do_bdist_wheel = True

    # mfem_source = './external/mfem'
    bglb.ext_prefix = os.path.join(bglb.prefix, 'mfem', 'external')
    print("ext_prefix(bdist)", bglb.ext_prefix)
    bglb.hypre_prefix = bglb.ext_prefix
    bglb.metis_prefix = bglb.ext_prefix

    bglb.mfem_prefix = bglb.ext_prefix
    bglb.mfems_prefix = os.path.join(bglb.ext_prefix, 'ser')
    bglb.mfemp_prefix = os.path.join(bglb.ext_prefix, 'par')

    bglb.mfem_build_miniapps = False


class Install(_install):
    '''
    called when pyton setup.py install
    '''
    user_options = _install.user_options + [
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
         'libmfem.so must exits under <mfems-prefix>/lib. ',
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
        ('git-sshclone', None, 'Use SSH for git clone',
         'try if default git clone using https fails (need Github account and setting for SSH)'),
        ('swig', None, 'Run Swig and exit'),
        ('skip-swig', None,
         'Skip running swig (used when wrapper is generated for the MFEM C++ library to be used'),
        ('ext-only', None, 'Build metis, hypre, mfem(C++) only'),
        ('skip-ext', None, 'Skip building metis, hypre, mfem(C++) only'),
        ('build-only', None, 'Skip final install stage to prefix'),
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

    def initialize_options(self):
        _install.initialize_options(self)
        self.swig = False
        self.skip_swig = False
        self.ext_only = False

        self.git_sshclone = False
        self.skip_ext = False
        self.with_parallel = False
        self.build_only = False
        self.no_serial = False
        self.mfem_prefix = ''
        self.mfems_prefix = ''
        self.mfemp_prefix = ''
        self.mfem_source = bglb.mfem_source
        self.mfem_branch = ''
        self.mfem_debug = False
        self.mfem_build_miniapps = False
        self.metis_prefix = ''
        self.hypre_prefix = ''

        self.with_cuda = False
        self.with_cuda_hypre = False
        self.cuda_arch = None
        self.with_metis64 = False

        self.with_pumi = False
        self.pumi_prefix = ''

        self.with_strumpack = False
        self.strumpack_prefix = ''

        self.with_suitesparse = False
        self.suitesparse_prefix = ''

        self.with_lapack = False
        self.blas_libraries = ""
        self.lapack_libraries = ""

        self.with_libceed = False
        self.libceed_prefix = ''
        self.libceed_only = False

        self.with_gslib = False
        self.gslib_prefix = ''
        self.gslib_only = False

        self.CC = ''
        self.CXX = ''
        self.MPICC = ''
        self.MPICXX = ''
        self.vv = False

        self.unverifiedSSL = False

    def finalize_options(self):
        if (bool(self.ext_only) and bool(self.skip_ext)):
            assert False, "skip-ext and ext-only can not use together"

        given_prefix = True
        if (self.prefix == '' or
                self.prefix is None):
            given_prefix = False
        else:
            self.prefix = os.path.expanduser(self.prefix)
            prefix = self.prefix

        bglb.verbose = bool(self.vv)
        if given_prefix:
            # global ext_prefix
            self.prefix = abspath(bglb.prefix)
            # ext_prefix = abspath(prefix)
        else:
            if '--user' in sys.argv:
                path = site.getusersitepackages()
                if not os.path.exists(path):
                    try:
                        print("attempting to make a --user directory", path)
                        os.makedirs(path)
                    except BaseException:
                        pass
                if os.path.exists(path):
                    path = os.path.dirname(path)
                    path = os.path.dirname(path)
                    path = os.path.dirname(path)
                else:
                    assert False, "no installation path found"
                self.prefix = path
            else:
                self.prefix = sys.prefix

        self.user = 0
        _install.finalize_options(self)

        bglb.use_unverifed_SSL = self.unverifiedSSL

        if bglb.verbose:
            print("prefix is :", self.prefix)

    def run(self):
        if not bglb.is_configured:
            configure_install(self)
            print_config()

        if bglb.swig_only:
            self.run_command("build")
        else:
            _install.run(self)


class BuildPy(_build_py):
    '''
    Called when python setup.py build_py
    '''
    user_options = _build_py.user_options

    def initialize_options(self):
        _build_py.initialize_options(self)

    def finalize_options(self):
        _build_py.finalize_options(self)

    def run(self):
        if not bglb.swig_only:
            if bglb.build_metis:
                if use_metis_gklib:
                    gitclone('gklib', use_sha=True)
                    gitclone('metis', use_sha=True)
                    make_metis(use_int64=metis_64, use_real64=metis_64)
                else:
                    download('metis')
                    make_metis(use_int64=metis_64, use_real64=metis_64)

            if bglb.build_hypre:
                download('hypre')
                cmake_make_hypre()
            if bglb.build_libceed:
                download('libceed')
                make_libceed()
            if bglb.build_gslib:
                download('gslib')
                make_gslib(serial=True)
                if bglb.build_hypre:
                    make_gslib()

            mfem_downloaded = False

            if bglb.build_mfem:
                gitclone('mfem', use_sha=True) if bglb.mfem_branch is None else gitclone(
                    'mfem', branch=bglb.mfem_branch)
                mfem_downloaded = True
                cmake_make_mfem(serial=True)

            if bglb.build_mfemp:
                if not mfem_downloaded:
                    gitclone('mfem', use_sha=True) if bglb.mfem_branch is None else gitclone(
                        'mfem', branch=bglb.mfem_branch)
                cmake_make_mfem(serial=False)

        if bglb.clean_swig:
            clean_wrapper()
        if bglb.run_swig:
            generate_wrapper()
            if bglb.swig_only:
                return

        if bglb.build_serial:
            make_mfem_wrapper(serial=True)
        if bglb.build_parallel:
            make_mfem_wrapper(serial=False)

        if not bglb.skip_install:
            _build_py.run(self)
        else:
            sys.exit()

class InstallLib(_install_lib):
    def finalize_options(self):
        _install_lib.finalize_options(self)
        src_cmd_obj = self.distribution.get_command_obj('install')
        src_cmd_obj.ensure_finalized()
        self.install_dir = src_cmd_obj.install_platlib


class InstallEggInfo(_install_egg_info):
    def run(self):
        if not bglb.dry_run:
            _install_egg_info.run(self)
        else:
            print("skipping regular install_egg_info")


class InstallScripts(_install_scripts):
    def run(self):
        if not bglb.dry_run:
            _install_scripts.run(self)
        else:
            print("skipping regular install_scripts")


class Clean(_clean):
    '''
    Called when python setup.py clean
    '''
    user_options = _clean.user_options + [
        ('ext', None, 'clean exteranal dependencies)'),
        ('mfem', None, 'clean mfem'),
        ('metis', None, 'clean metis'),
        ('hypre', None, 'clean hypre'),
        ('swig', None, 'clean swig'),
        ('all-exts', None, 'delete all externals'),
    ]

    def initialize_options(self):
        _clean.initialize_options(self)
        self.ext = False
        self.mfem = False
        self.hypre = False
        self.metis = False
        self.swig = False
        self.all_exts = False

    def run(self):
        bglb.dry_run = self.dry_run
        bglb.verbose = bool(self.verbose)

        os.chdir(extdir)

        make_command = find_command('make')

        if self.ext or self.mfem:
            path = os.path.join(extdir, 'mfem', 'cmbuild_par')
            if os.path.exists(path):
                shutil.rmtree(path)
            path = os.path.join(extdir, 'mfem', 'cmbuild_ser')
            if os.path.exists(path):
                shutil.rmtree(path)
        if self.ext or self.hypre:
            path = os.path.join(extdir, 'hypre', 'cmbuild')
            if os.path.exists(path):
                shutil.rmtree(path)
        if self.ext or self.metis:
            path = os.path.join(extdir, 'metis')
            if os.path.exists(path):
                os, chdir(path)
                command = ['make', 'clean']
                subprocess.check_call(command)
        if self.all_exts or self.all:
            for xxx in ('metis', 'hypre', 'mfem', 'gslib', 'gklib', 'libceed'):
                path = os.path.join(extdir, xxx)
                if os.path.exists(path):
                    shutil.rmtree(path)

        if self.swig or self.all:
            clean_wrapper()

        clean_so(all=self.all)

        os.chdir(rootdir)
        _clean.run(self)


if __name__ == '__main__':
    cmdclass = {'build_py': BuildPy,
                'install': Install,
                'install_lib': InstallLib,
                'install_egg_info': InstallEggInfo,
                'install_scripts': InstallScripts,
                'clean': Clean}
    if haveWheel:
        cmdclass['bdist_wheel'] = BdistWheel

    setup(
        cmdclass=cmdclass,
    )
