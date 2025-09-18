# ----------------------------------------------------------------------------------------
# Routines for PyMFEM Wrapper Generation/Compile
# ----------------------------------------------------------------------------------------
import sys
import os
import re
import subprocess

__all__ = ["cmake_make_mfem"]

from build_consts import *
from build_utils import *

import build_globals as bglb


def cmake_make_mfem(serial=True):
    '''
    build MFEM
    '''
    cmbuild = 'cmbuild_ser' if serial else 'cmbuild_par'
    path = os.path.join(extdir, 'mfem', cmbuild)
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)

    ldflags = os.getenv('LDFLAGS') if os.getenv('LDFLAGS') is not None else ''
    metisflags = ''
    hypreflags = ''

    rpaths = []
    if sys.platform in ("linux", "linux2"):
        rpaths_origin = "$ORIGIN"
    elif sys.platform == "darwin":
        rpaths_origin = "@loader_path"

    def add_rpath(p, dest):
        p = os.path.join(rpaths_origin, os.path.relpath(p, dest))
        
        if not p in rpaths:
            rpaths.append(p)

    cmake_opts = {'DBUILD_SHARED_LIBS': '1',
                  'DMFEM_ENABLE_EXAMPLES': '1',
                  'DMFEM_ENABLE_MINIAPPS': '0',
                  'DCMAKE_SHARED_LINKER_FLAGS': ldflags,
                  'DMFEM_USE_ZLIB': '1',
                  'DCMAKE_CXX_FLAGS': bglb.cxxstd_flag,
                  'DCMAKE_BUILD_WITH_INSTALL_RPATH': '1'}

    if sys.platform == 'darwin':
        cmake_opts["DCMAKE_MACOSX_RPATH"] = 'YES'
        cmake_opts["DCMAKE_INSTALL_NAME_DIR"] = '@rpath'        

    if bglb.mfem_debug:
        cmake_opts['DMFEM_DEBUG'] = 'YES'

    if bglb.mfem_build_miniapps:
        cmake_opts['DMFEM_ENABLE_MINIAPPS'] = '1'

    if bglb.verbose:
        cmake_opts['DCMAKE_VERBOSE_MAKEFILE'] = '1'

    if serial:
        ex_loc = os.path.join(bglb.mfems_prefix, "examples")        
        cmake_opts['DCMAKE_CXX_COMPILER'] = bglb.cxx_command
        cmake_opts['DMFEM_USE_EXCEPTIONS'] = '1'
        cmake_opts['DCMAKE_INSTALL_PREFIX'] = bglb.mfems_prefix

        add_rpath(os.path.join(bglb.mfems_prefix, 'lib'), ex_loc)
        if bglb.enable_suitesparse:
            enable_metis = True
        else:
            enable_metis = False
        #assert False, rpaths            
    else:
        ex_loc = os.path.join(bglb.mfemp_prefix, "examples")
        cmake_opts['DCMAKE_CXX_COMPILER'] = bglb.mpicxx_command
        cmake_opts['DMFEM_USE_EXCEPTIONS'] = '0'
        cmake_opts['DCMAKE_INSTALL_PREFIX'] = bglb.mfemp_prefix
        cmake_opts['DMFEM_USE_MPI'] = '1'
        cmake_opts['DHYPRE_DIR'] = bglb.hypre_prefix
        cmake_opts['DHYPRE_INCLUDE_DIRS'] = os.path.join(
            bglb.hypre_prefix, "include")

        add_rpath(os.path.join(bglb.mfemp_prefix, 'lib'), ex_loc)

        hyprelibpath = os.path.dirname(
            find_libpath_from_prefix(
                'HYPRE', bglb.hypre_prefix))

        add_rpath(hyprelibpath, ex_loc)

        hypreflags = "-L" + hyprelibpath + " -lHYPRE "

        if bglb.enable_strumpack:
            cmake_opts['DMFEM_USE_STRUMPACK'] = '1'
            cmake_opts['DSTRUMPACK_DIR'] = bglb.strumpack_prefix
            libpath = os.path.dirname(
                find_libpath_from_prefix("STRUMPACK", bglb.strumpack_prefix))
            add_rpath(libpath, ex_loc)
        if bglb.enable_pumi:
            cmake_opts['DMFEM_USE_PUMI'] = '1'
            cmake_opts['DPUMI_DIR'] = bglb.pumi_prefix
            libpath = os.path.dirname(
                find_libpath_from_prefix("pumi", bglb.strumpack_prefix))
            add_rpath(libpath, ex_loc)
        enable_metis = True

    if enable_metis:
        cmake_opts['DMFEM_USE_METIS_5'] = '1'
        cmake_opts['DMETIS_DIR'] = bglb.metis_prefix
        cmake_opts['DMETIS_INCLUDE_DIRS'] = os.path.join(
            bglb.metis_prefix, "include")
        metislibpath = os.path.dirname(
            find_libpath_from_prefix(
                'metis', bglb.metis_prefix))
        add_rpath(metislibpath, ex_loc)

        if bglb.use_metis_gklib:
            metisflags = "-L" + metislibpath + " -lmetis -lGKlib "
        else:
            metisflags = "-L" + metislibpath + " -lmetis "

    if ldflags != '':
        cmake_opts['DCMAKE_SHARED_LINKER_FLAGS'] = ldflags
        cmake_opts['DCMAKE_EXE_LINKER_FLAGS'] = ldflags

    if metisflags != '':
        cmake_opts['DMETIS_LIBRARIES'] = metisflags
    if hypreflags != '':
        cmake_opts['DHYPRE_LIBRARIES'] = hypreflags

    if bglb.enable_cuda:
        cmake_opts['DMFEM_USE_CUDA'] = '1'
        if bglb.cuda_arch != '':
            cmake_opts['DCMAKE_CUDA_ARCHITECTURES'] = bglb.cuda_arch

    if bglb.enable_libceed:
        cmake_opts['DMFEM_USE_CEED'] = '1'
        cmake_opts['DCEED_DIR'] = bglb.libceed_prefix
        libpath = os.path.dirname(
            find_libpath_from_prefix("ceed", bglb.libceed_prefix))
        add_rpath(libpath, ex_loc)

    if bglb.enable_gslib:
        if serial:
            cmake_opts['DMFEM_USE_GSLIB'] = '1'
            cmake_opts['DGSLIB_DIR'] = bglb.gslibs_prefix
        else:
            cmake_opts['DMFEM_USE_GSLIB'] = '1'
            cmake_opts['DGSLIB_DIR'] = bglb.gslibp_prefix

    if bglb.enable_suitesparse:
        cmake_opts['DMFEM_USE_SUITESPARSE'] = '1'
        if bglb.suitesparse_prefix != '':
            cmake_opts['DSuiteSparse_DIR'] = bglb.suitesparse_prefix

    if bglb.enable_lapack:
        cmake_opts['DMFEM_USE_LAPACK'] = '1'
    if bglb.blas_libraries != "":
        cmake_opts['DBLAS_LIBRARIES'] = bglb.blas_libraries
    if bglb.lapack_libraries != "":
        cmake_opts['DLAPACK_LIBRARIES'] = bglb.lapack_libraries

    cmake_opts['DCMAKE_INSTALL_RPATH'] = ";".join(rpaths)

    pwd = chdir(path)
    cmake('..', **cmake_opts)

    txt = 'serial' if serial else 'parallel'

    make('mfem_' + txt)
    make_install('mfem_' + txt)

    from shutil import copytree, rmtree

    print("copying mesh data for testing", "../data",
          cmake_opts['DCMAKE_INSTALL_PREFIX'])
    path = os.path.join(cmake_opts['DCMAKE_INSTALL_PREFIX'], "data")
    if os.path.exists(path):
        rmtree(path)
    copytree("../data", path)

    os.chdir(pwd)
