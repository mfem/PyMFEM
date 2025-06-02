# ----------------------------------------------------------------------------------------
# Routines for libceed
# ----------------------------------------------------------------------------------------
import sys
import os
import re
import subprocess

from build_utils import (
    read_mfem_tplflags, abspath, external_install_prefix,
    make_call, chdir, remove_files, download, gitclone,
    record_mfem_sha, get_numpy_inc, get_mpi4py_inc, find_libpath_from_prefix,]

import build_globals as bglb
from build_consts as import *

__all__ = ["make_libceed"]


def make_libceed(serial=False):
    if bglb.verbose:
        print("Building libceed")

    path = os.path.join(extdir, 'libceed')
    if not os.path.exists(path):
        assert False, "libceed is not downloaded"

    pwd = chdir(path)
    try:
        make_call(['make', 'clean'])
    except:
        pass

    if enable_cuda:
        command = ['make', 'configure', 'CUDA_DIR='+cuda_prefix]
        make_call(command)

    make('libceed')
    make_install('libceed', prefix=libceed_prefix)
    os.chdir(pwd)
