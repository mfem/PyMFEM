# ----------------------------------------------------------------------------------------
# Routines for libceed
# ----------------------------------------------------------------------------------------
import sys
import os
import re
import subprocess

import build_globals as bglb
from build_consts import *
from build_utils import *

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

    if bglb.enable_cuda:
        command = ['make', 'configure', 'CUDA_DIR='+bglb.cuda_prefix]
        make_call(command)

    make('libceed')
    make_install('libceed', prefix=bglb.libceed_prefix)
    os.chdir(pwd)
