# ----------------------------------------------------------------------------------------
# Routines for metis
# ----------------------------------------------------------------------------------------
import sys
import os
import re
import subprocess
from sys import platform

__all__ = ["make_metis_gklib", "make_metis"]

from build_utils import *
from build_consts import *

import build_globals as bglb


def make_metis_gklib(use_int64=False, use_real64=False):
    '''
    build GKlib/metis
    '''

    '''
    build/install GKlib
    '''
    if bglb.verbose:
        print("Building gklib")

    path = os.path.join(extdir, 'gklib')
    if not bglb.dry_run and not os.path.exists(path):
        assert False, "gklib is not downloaded"

    path = os.path.join(path, 'cmbuild')
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)
    pwd = chdir(path)

    cmake_opts = {'DBUILD_SHARED_LIBS': '1',
                  'DCMAKE_INSTALL_PREFIX': metis_prefix}
    if verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    cmake('..', **cmake_opts)
    make('gklib')
    make_install('gklib')
    os.chdir(pwd)

    '''
    build/install metis
    '''
    path = os.path.join(extdir, 'metis')
    if not bglb.dry_run and not os.path.exists(path):
        assert False, "metis is not downloaded"
    elif not os.path.exists(path):
        os.makedirs(path)
        os.makedirs(os.path.join(path, 'build'))

    pwd = chdir(path)

    gklibpath = os.path.dirname(find_libpath_from_prefix(
        'GKlib', bglb.metis_prefix))

    options = ['gklib_path='+bglb.metis_prefix]
    if use_int64:
        options.append('i64=1')

    if use_real64:
        options.append('r64=1')

    command = ['make', 'config', 'shared=1'] + options
    command = command + ['prefix=' + bglb.metis_prefix, 'cc=' + bglb.cc_command]
    make_call(command)

    chdir('build')
    cmake_opts = {'DGKLIB_PATH': bglb.metis_prefix,
                  'DSHARED': '1',
                  'DCMAKE_C_COMPILER': bglb.cc_command,
                  'DCMAKE_C_STANDARD_LIBRARIES': '-lGKlib',
                  'DCMAKE_INSTALL_RPATH': "@loader_path",
                  'DCMAKE_BUILD_WITH_INSTALL_RPATH': '1',
                  'DCMAKE_INSTALL_PREFIX': bglb.metis_prefix}
    if bglb.verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    cmake('..', **cmake_opts)
    chdir(path)
    make('metis')
    make_install('metis')

    if platform == "darwin":
        command = ['install_name_tool',
                   '-id',
                   os.path.join("@rpath", 'libGKlib.dylib'),
#                   os.path.join(bglb.metis_prefix, 'lib', 'libGKlib.dylib'),
                   os.path.join(bglb.metis_prefix, 'lib', 'libGKlib.dylib'), ]
        make_call(command)
        command = ['install_name_tool',
                   '-id', 
                   os.path.join("@rpath", 'libmetis.dylib'),
                   os.path.join(bglb.metis_prefix, 'lib', 'libmetis.dylib'), ]
        make_call(command)
    os.chdir(pwd)


def make_metis(use_int64=False, use_real64=False):
    '''
    build metis
    '''
    if bglb.verbose:
        print("Building metis")

    path = os.path.join(extdir, 'metis')
    if not os.path.exists(path):
        assert False, "metis is not downloaded"

    pwd = chdir(path)

    if use_int64:
        pattern_int = "#define IDXTYPEWIDTH 32"
        replace_int = "#define IDXTYPEWIDTH 64"
    else:
        pattern_int = "#define IDXTYPEWIDTH 64"
        replace_int = "#define IDXTYPEWIDTH 32"
    with open("include/metis.h", "r") as metis_header_fid:
        metis_header_lines = metis_header_fid.readlines()
    with open("include/metis.h", "w") as metis_header_fid:
        for line in metis_header_lines:
            metis_header_fid.write(re.sub(pattern_int, replace_int, line))

    if use_real64:
        pattern_real = "#define REALTYPEWIDTH 32"
        replace_real = "#define REALTYPEWIDTH 64"
    else:
        pattern_real = "#define REALTYPEWIDTH 64"
        replace_real = "#define REALTYPEWIDTH 32"
    with open("include/metis.h", "r") as metis_header_fid:
        metis_header_lines = metis_header_fid.readlines()
    with open("include/metis.h", "w") as metis_header_fid:
        for line in metis_header_lines:
            metis_header_fid.write(re.sub(pattern_real, replace_real, line))

    command = ['make', 'config', 'shared=1',
               'prefix=' + bglb.metis_prefix,
               'cc=' + bglb.cc_command]
    make_call(command, env={'CMAKE_POLICY_VERSION_MINIMUM': '3.5'})
    make('metis')
    make_install('metis')

    if platform == "darwin":
        command = ['install_name_tool',
                   '-id',
                   os.path.join('@rpath', 'libmetis.dylib'),
                   os.path.join(bglb.metis_prefix, 'lib', 'libmetis.dylib'), ]
        make_call(command)
    os.chdir(pwd)
