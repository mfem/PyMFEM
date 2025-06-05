# ----------------------------------------------------------------------------------------
# Routines for hypre
# ----------------------------------------------------------------------------------------

import sys
import os
import re
import subprocess

__all__ = ["cmake_make_hypre"]

from build_utils import *
from build_consts import *

import build_globals as bglb


def cmake_make_hypre():
    '''
    build hypre
    '''
    if bglb.verbose:
        print("Building hypre")

    cmbuild = 'cmbuild'
    path = os.path.join(extdir, 'hypre', 'src', cmbuild)
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)

    pwd = chdir(path)

    cmake_opts = {'DBUILD_SHARED_LIBS': '1',
                  'DHYPRE_INSTALL_PREFIX': bglb.hypre_prefix,
                  'DHYPRE_ENABLE_SHARED': '1',
                  'DCMAKE_C_FLAGS': '-fPIC',
                  'DCMAKE_INSTALL_PREFIX': bglb.hypre_prefix,
                  'DCMAKE_INSTALL_NAME_DIR': "@rpath", }
    if bglb.verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    if bglb.enable_cuda and bglb.enable_cuda_hypre:
        # in this case, settitng CMAKE_C_COMPILER
        # causes "mpi.h" not found error. For now, letting CMAKE
        # to find MPI
        cmake_opts['DHYPRE_WITH_CUDA'] = '1'
        if bglb.cuda_arch != '':
            cmake_opts['DCMAKE_CUDA_ARCHITECTURES'] = bglb.cuda_arch
    else:
        cmake_opts['DCMAKE_C_COMPILER'] = bglb.mpicc_command

    cmake('..', **cmake_opts)

    make('hypre')
    make_install('hypre')

    os.chdir(pwd)
