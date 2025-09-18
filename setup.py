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
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
from distutils.command.clean import clean as _clean

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_build_system"))

import build_globals as bglb
from build_mfem import *
from build_metis import *
from build_hypre import *
from build_pymfem import *
from build_libceed import *
from build_gslib import *
from build_config import *
from build_consts import *
from build_utils import *

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------

class Install(_install):
    '''
    called when pyton setup.py install
    '''
    user_options = _install.user_options + cmd_options

    def initialize_options(self):
        _install.initialize_options(self)
        initialize_cmd_options(self)

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
                        print("!!!!! attempting to make a --user directory", path)
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
            print("!!!!! prefix is :", self.prefix)

    def run(self):
        if not bglb.is_configured:
            if bglb.verbose:
                 print('!!!!! Running config (install)')
            bglb.prefix = abspath(self.prefix)
            configure_build(self)
            print_config()

        if bglb.swig_only and not bglb.build_py_done:
            #  comes here if python setup.py install is used
            self.run_command("build")
        else:
            _install.run(self)


class BdistWheel(_bdist_wheel):
    '''
    Wheel build performs SWIG + Serial in Default.
    --skip-build option skip building entirely.
    '''
    user_options = _bdist_wheel.user_options + cmd_options

    def initialize_options(self):
        _bdist_wheel.initialize_options(self)
        initialize_cmd_options(self)

    def finalize_options(self):
        def _has_ext_modules():
            return True
        self.distribution.has_ext_modules = _has_ext_modules
        _bdist_wheel.finalize_options(self)

    def run(self):
        import build_globals as bglb

        if not bglb.is_configured:
            if bglb.verbose:
                print('!!!!! Running config (bdist wheel)')
            bglb.prefix = abspath(self.bdist_dir)
            bglb.ext_prefix = os.path.join(bglb.prefix, 'mfem', 'external')
            bglb.bdist_wheel_dir = abspath(self.bdist_dir)
            bglb.do_bdist_wheel = True

            configure_build(self)
            clean_dist_info(bglb.prefix)
            if bglb.keep_temp:
                self.keep_temp = True
            print_config()

        self.run_command("build")

        _bdist_wheel.run(self)


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
        bglb.build_py_done = True

        if not bglb.swig_only:
            if bglb.build_metis:
                if bglb.use_metis_gklib:
                    gitclone('gklib', use_sha=True)
                    gitclone('metis', use_sha=True)
                    make_metis(use_int64=bglb.metis_64,
                               use_real64=bglb.metis_64)
                else:
                    download('metis')
                    make_metis(use_int64=bglb.metis_64,
                               use_real64=bglb.metis_64)

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
            generate_wrapper(bglb.run_swig_parallel)

        if bglb.build_serial:
            make_mfem_wrapper(serial=True)
        if bglb.build_parallel:
            make_mfem_wrapper(serial=False)

        _build_py.run(self)

class InstallLib(_install_lib):
    def finalize_options(self):
        _install_lib.finalize_options(self)
        src_cmd_obj = self.distribution.get_command_obj('install')
        src_cmd_obj.ensure_finalized()
        self.install_dir = src_cmd_obj.install_platlib

    def run(self):
        if not bglb.dry_run:
            _install_lib.run(self)
        else:
            print("skipping regular install_lib")


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
                'clean': Clean,
                'bdist_wheel': BdistWheel}

    setup(
        cmdclass=cmdclass,
    )
