# ----------------------------------------------------------------------------------------
# Routines for PyMFEM Wrapper Generation/Compile
# ----------------------------------------------------------------------------------------
import sys
import os
import re
import subprocess

__all__ = ["write_setup_local", "generate_wrapper",
           "clean_wrapper", "make_mfem_wrapper"]

from build_utils import *
from build_consts import *
import build_globals as bglb


def write_setup_local():
    '''
    create setup_local.py. parameters written here will be read
    by setup.py in mfem._ser and mfem._par
    '''
    mfemser = bglb.mfems_prefix
    mfempar = bglb.mfemp_prefix

    hyprelibpath = os.path.dirname(
        find_libpath_from_prefix('HYPRE', bglb.hypre_prefix))
    metislibpath = os.path.dirname(
        find_libpath_from_prefix('metis', bglb.metis_prefix))

    mfems_tpl = read_mfem_tplflags(bglb.mfems_prefix)
    mfemp_tpl = read_mfem_tplflags(
        bglb.mfemp_prefix) if bglb.build_parallel else ''

    print(mfems_tpl, mfemp_tpl)

    params = {'cxx_ser': bglb.cxx_command,
              'cc_ser': bglb.cc_command,
              'cxx_par': bglb.mpicxx_command,
              'cc_par': bglb.mpicc_command,
              'whole_archive': '--whole-archive',
              'no_whole_archive': '--no-whole-archive',
              'nocompactunwind': '',
              'swigflag': '-Wall -c++ -python -fastproxy -olddefs -keyword',
              'hypreinc': os.path.join(bglb.hypre_prefix, 'include'),
              'hyprelib': hyprelibpath,
              'metisinc': os.path.join(bglb.metis_prefix, 'include'),
              'metis5lib': metislibpath,
              'numpyinc': get_numpy_inc(),
              'mpi4pyinc': '',
              'mpiinc':bglb.mpiinc,
              'mfem_outside': '1' if  bglb.mfem_outside else '0',
              'mfembuilddir': os.path.join(mfempar, 'include'),
              'mfemincdir': os.path.join(mfempar, 'include', 'mfem'),
              'mfemlnkdir': os.path.join(mfempar, 'lib'),
              'mfemserbuilddir': os.path.join(mfemser, 'include'),
              'mfemserincdir': os.path.join(mfemser, 'include', 'mfem'),
              'mfemserlnkdir': os.path.join(mfemser, 'lib'),
              'mfemsrcdir': os.path.join(bglb.mfem_source),
              'mfemstpl': mfems_tpl,
              'mfemptpl': mfemp_tpl,
              'add_pumi': '',
              'add_strumpack': '',
              'add_cuda': '',
              'add_libceed': '',
              'add_suitesparse': '',
              'add_gslib': '',
              'add_gslibp': '',
              'add_gslibs': '',
              'libceedinc': os.path.join(bglb.libceed_prefix, 'include'),
              'gslibsinc': os.path.join(bglb.gslibs_prefix, 'include'),
              'gslibpinc': os.path.join(bglb.gslibp_prefix, 'include'),
              'cxxstdflag': bglb.cxxstd_flag,
              'build_mfem': '1' if bglb.build_mfem else '0',
              'bdist_wheel_dir': bglb.bdist_wheel_dir,
              }

    if bglb.build_parallel:
        params['mpi4pyinc'] = get_mpi4py_inc()

    def add_extra(xxx, inc_sub=None):
        params['add_' + xxx] = '1'
        ex_prefix = getattr(bglb, xxx + '_prefix')
        if inc_sub is None:
            params[xxx +
                   'inc'] = os.path.join(ex_prefix, 'include')
        else:
            params[xxx +
                   'inc'] = os.path.join(ex_prefix, 'include', inc_sub)

        params[xxx + 'lib'] = os.path.join(ex_prefix, 'lib')

    if bglb.enable_pumi:
        add_extra('pumi')
    if bglb.enable_strumpack:
        add_extra('strumpack')
    if bglb.enable_cuda:
        add_extra('cuda')
    if bglb.enable_libceed:
        add_extra('libceed')
    if bglb.enable_suitesparse:
        add_extra('suitesparse', inc_sub='suitesparse')
    if bglb.enable_gslib:
        add_extra('gslibs')
    if bglb.enable_gslib:
        add_extra('gslibp')

    pwd = chdir(rootdir)

    fid = open('setup_local.py', 'w')
    fid.write("#  setup_local.py \n")
    fid.write("#  generated from setup.py\n")
    fid.write("#  do not edit this directly\n")

    for key, value in params.items():
        text = key.lower() + ' = "' + value + '"'
        fid.write(text + "\n")
    fid.close()

    os.chdir(pwd)


def generate_wrapper(do_parallel):
    '''
    run swig.
    '''
    # this should work as far as we are in the same directory ?
    from multiprocessing import Pool, cpu_count
    import build_globals as bglb

    if bglb.dry_run or bglb.verbose:
        print("generating SWIG wrapper")
        print("using MFEM source", os.path.abspath(bglb.mfem_source))
    if not os.path.exists(os.path.abspath(bglb.mfem_source)):
        assert False, "MFEM source directory. Use --mfem-source=<path>"

    def ifiles():
        ifiles = os.listdir()
        ifiles = [x for x in ifiles if x.endswith('.i')]
        ifiles = [x for x in ifiles if not x.startswith('#')]
        ifiles = [x for x in ifiles if not x.startswith('.')]
        return ifiles

    def check_new(ifile):
        wfile = ifile[:-2] + '_wrap.cxx'
        if not os.path.exists(wfile):
            return True
        return os.path.getmtime(ifile) > os.path.getmtime(wfile)

    def update_integrator_exts():
        pwd = chdir(os.path.join(rootdir, 'mfem', 'common'))
        command1 = [sys.executable, "generate_lininteg_ext.py"]
        command2 = [sys.executable, "generate_bilininteg_ext.py"]
        make_call(command1)
        make_call(command2)
        os.chdir(pwd)

    def update_header_exists(mfem_source):
        print("updating the list of existing headers")
        list_of_headers = []
        L = len(mfem_source.split(os.sep))
        for (dirpath, dirnames, filenames) in os.walk(mfem_source):
            for filename in filenames:
                if filename.endswith('.hpp'):
                    dirs = dirpath.split(os.sep)[L:]
                    dirs.append(filename[:-4])
                    tmp = '_'.join(dirs)
                    xx = re.split('_|-', tmp)
                    new_name = 'FILE_EXISTS_'+'_'.join([x.upper() for x in xx])
                    if new_name not in list_of_headers:
                        list_of_headers.append(new_name)

        pwd = chdir(os.path.join(rootdir, 'mfem', 'common'))
        fid = open('existing_mfem_headers.i', 'w')
        for x in list_of_headers:
            fid.write("#define " + x + "\n")
        fid.close()
        os.chdir(pwd)

    mfemser = bglb.mfems_prefix
    mfempar = bglb.mfemp_prefix

    update_header_exists(bglb.mfem_source)

    swigflag = '-Wall -c++ -python -fastproxy -olddefs -keyword'.split(' ')

    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))

    serflag = ['-I' + os.path.join(mfemser, 'include'),
               '-I' + os.path.join(mfemser, 'include', 'mfem'),
               '-I' + os.path.abspath(bglb.mfem_source)]
    if bglb.enable_suitesparse:
        serflag.append('-I' + os.path.join(bglb.suitesparse_prefix,
                                           'include', 'suitesparse'))

    for filename in ['lininteg.i', 'bilininteg.i']:
        command = [swig_command] + swigflag + serflag + [filename]
        make_call(command)
    update_integrator_exts()

    commands = []
    for filename in ifiles():
        if not check_new(filename):
            continue
        command = [swig_command] + swigflag + serflag + [filename]
        commands.append(command)

    mp_pool = Pool(max((cpu_count() - 1, 1)))
    with mp_pool:
        mp_pool.map(subprocess.run, commands)

    if not do_parallel:
        os.chdir(pwd)
        return

    chdir(os.path.join(rootdir, 'mfem', '_par'))

    parflag = ['-I' + os.path.join(mfempar, 'include'),
               '-I' + os.path.join(mfempar, 'include', 'mfem'),
               '-I' + os.path.abspath(bglb.mfem_source),
               '-I' + os.path.join(bglb.hypre_prefix, 'include'),
               '-I' + os.path.join(bglb.metis_prefix, 'include'),
               '-I' + get_mpi4py_inc()]

    if bglb.enable_pumi:
        parflag.append('-I' + os.path.join(bglb.pumi_prefix, 'include'))
    if bglb.enable_strumpack:
        parflag.append('-I' + os.path.join(bglb.strumpack_prefix, 'include'))
    if bglb.enable_suitesparse:
        parflag.append('-I' + os.path.join(bglb.suitesparse_prefix,
                                           'include', 'suitesparse'))

    commands = []
    for filename in ifiles():
        if filename == 'strumpack.i' and not bglb.enable_strumpack:
            continue
        if not check_new(filename):
            continue
        command = [swig_command] + swigflag + parflag + [filename]
        commands.append(command)

    mp_pool = Pool(max((cpu_count() - 1, 1)))
    with mp_pool:
        mp_pool.map(subprocess.run, commands)

    os.chdir(pwd)


def clean_wrapper():
    from pathlib import Path

    # serial
    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.cxx')]

    print(os.getcwd(), wfiles)
    remove_files(wfiles)

    wfiles = [x for x in os.listdir() if x.endswith('_wrap.h')]
    remove_files(wfiles)

    wfiles = [x for x in os.listdir() if x.endswith('.py')]
    wfiles.remove("__init__.py")
    wfiles.remove("setup.py")
    wfiles.remove("tmop_modules.py")
    remove_files(wfiles)

    ifiles = [x for x in os.listdir() if x.endswith('.i')]
    for x in ifiles:
        Path(x).touch()

    # parallel
    chdir(os.path.join(rootdir, 'mfem', '_par'))
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.cxx')]

    remove_files(wfiles)
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.h')]
    remove_files(wfiles)

    wfiles = [x for x in os.listdir() if x.endswith('.py')]
    wfiles.remove("__init__.py")
    wfiles.remove("setup.py")
    wfiles.remove("tmop_modules.py")
    remove_files(wfiles)

    ifiles = [x for x in os.listdir() if x.endswith('.i')]
    for x in ifiles:
        Path(x).touch()

    chdir(pwd)


def make_mfem_wrapper(serial=True):
    '''
    compile PyMFEM wrapper code
    '''
    from multiprocessing import cpu_count
    import build_globals as bglb

    if bglb.dry_run or bglb.verbose:
        print("compiling wrapper code, serial=" + str(serial))
    if not os.path.exists(os.path.abspath(bglb.mfem_source)):
        assert False, "MFEM source directory. Use --mfem-source=<path>"

    record_mfem_sha(bglb.mfem_source)

    write_setup_local()

    if serial:
        pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))
    else:
        pwd = chdir(os.path.join(rootdir, 'mfem', '_par'))

    python = sys.executable
    command = [python, 'setup.py', 'build_ext', '--inplace', '--parallel',
               str(max((cpu_count() - 1, 1)))]
    make_call(command, force_verbose=True)

    os.chdir(pwd)
