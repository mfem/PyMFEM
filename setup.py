"""

  MFEM setuptools based setup module.
  

  python setup.py build  # build mfem and PyMFEM in serial
  python setup.py build --parallel # build metis/hypre/mfem and PyMFEM in parallel

  # choosing compiler
  python setup.py build --parallel --CC=xxx, --CXX=xxx, --MPICC=xxx, --MPICXX=xxx

  (plan)
  python setup.py build --cuda

  (note we will install every externals under <prefix>/mfem/)
     <prefix> /mfem/par  : mfem parallel
     <prefix> /mfem/ser  : mfem serial 
     <prefix> /mfem/lib  : mfem serial 
     <prefix> /mfem/lib  : mfem serial 
"""
import sys
import os
import urllib
import gzip
import re
import shutil
import subprocess
from subprocess import DEVNULL 

import multiprocessing

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.install import install as _install
try:
    from setuptools._distutils.command.clean import clean as _clean
except ImportError:
    from distutils.command.clean import clean as _clean

# To use a consistent encoding
from codecs import open

### constants
repos = {"mfem": "https://github.com/mfem/mfem/archive/v4.2.tar.gz",
         "metis": "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz",
         "hypre": "https://github.com/hypre-space/hypre/archive/v2.20.0.tar.gz",}

rootdir = os.path.abspath(os.path.dirname(__file__))
extdir = os.path.join(rootdir, 'external')

from sys import platform
if platform == "linux" or platform == "linux2":
    dylibext = '.so'
elif platform == "darwin":
    # OS X
    dylibext = '.dylib'
elif platform == "win32":
    # Windows...
    assert False, "Windows is not supported yet. Contribution is welcome"


### global variables
prefix = ''
verbose = False
build_parallel = False
build_serial = False
mfem_prefix = ''
metis_prefix = ''
hypre_prefix = ''

enable_pumi = False
pumi_prefix = ''

dry_run = False

cc_command = 'cc' if os.getenv("CC") is None else os.getenv("CC")
cxx_command = 'c++' if os.getenv("CC") is None else os.getenv("CXX")
mpicc_command = 'mpicc' if os.getenv("MPICC") is None else os.getenv("MPICC")
mpicxx_command = 'mpic++' if os.getenv("MPICXX") is None else os.getenv("MPICXX")
cxx11_flag = '-std=c++11'


### meta data
def version():
    VERSIONFILE = os.path.join('mfem', '__init__.py')
    initfile_lines = open(VERSIONFILE, 'rt').readlines()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    for line in initfile_lines:
        mo = re.search(VSRE, line, re.M)
        if mo:
            return mo.group(1)
    raise RuntimeError('Unable to find version string in %s.' % (VERSIONFILE,))

def long_description():
    with open(os.path.join(rootdir, 'README')) as f:
        return f.read()

keywords = """
scientific computing
finite element method
"""

platforms = """
Mac OS X
Linux
"""

metadata = {
    'name'             : 'mfem',
    'version'          : version(),
    'description'      : __doc__.strip(),
    'long_description' : long_description(),
    'url'              : 'http://mfem.org',
    'download_url'     : 'https://github.com/mfem',
    'classifiers'      :[#   3 - Alpha
                         #   4 - Beta
                         #   5 - Production/Stable
                         'Development Status :: 4 - Beta',
                         'Intended Audience :: Developers',
                         'Topic :: Scientific/Engineering :: Physics'
                         'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
                         'Programming Language :: Python :: 3.6',
                        ],
    
    'keywords'         : [k for k in keywords.split('\n')    if k],
    'platforms'        : [p for p in platforms.split('\n')   if p],
    'license'          : 'LGPL-2.1',    
    'author'           : 'MFEM developement team',
    'author_email'     : '',
    'maintainer'       : 'S. Shiraiwa',
    'maintainer_email' : 'shiraiwa@princeton.edu',
    }

##  utilities
def find_command(name):
    from shutil import which
    return which(name)

def make_call(command, target=''):
    '''
    call command
    '''
    if dry_run:
        print("calling ... " + " ".join(command))
        return
    kwargs = {}
    if not verbose:
        kwargs['stdout'] = DEVNULL
        kwargs['stderr'] = DEVNULL
    try:
        subprocess.check_call(command, **kwargs)
    except subprocess.CalledProcessError:
        if target == '':
            target = " ".join(command)
        print("Failed when calling command: " + target)
        raise

def chdir(path):
    '''
    change directory
    '''
    pwd = os.getcwd()
    os.chdir(path)
    if verbose:
        print("Moving to a directory : " + path)
    return pwd

def remove_files(files):
    for f in files:
        if verbose:
            print("Removing : " + f)
        if dry_run:
            continue
        os.remove(f)
    
def make(target):
    '''
    make : add -j option automatically
    '''
    command = ['make', '-j', str(max((multiprocessing.cpu_count()-1, 1)))]
    make_call(command, target=target) 

def make_install(target):
    '''
    make install
    '''
    command = ['make', 'install']
    make_call(command, target=target)

def download(xxx):
    '''
    download tar.gz from somewhere. xxx is name.
    url is given by repos above
    '''
    from urllib import request
    import tarfile

    if os.path.exists(os.path.join(extdir, xxx)):
        print("Download " + xxx + " skipped. Use clean --all-externals if needed")
        return
    url = repos[xxx]
    print("Downloading :", url)
    
    ftpstream = request.urlopen(url)
    targz = tarfile.open(fileobj=ftpstream, mode="r|gz")
    targz.extractall(path=extdir)
    os.rename(os.path.join(extdir, targz.getnames()[0].split('/')[0]),
              os.path.join(extdir, xxx))

def cmake(path, **kwargs):
    '''
    run cmake. must be called in the target directory
    '''
    command = ['cmake', path]
    for key, value in kwargs.items():
        command.append('-'+key+'='+value)
    make_call(command)

def cmake_make_hypre():
    '''
    build hypre
    '''
    cmbuild = 'cmbuild'
    path = os.path.join(extdir, 'hypre', 'src', cmbuild)
    if os.path.exists(path):
        print("working directory already exists!")
    else:
        os.makedirs(path)

    pwd = chdir(path)

    cmake_opts = {'DCMAKE_VERBOSE_MAKEFILE':'1',
                  'DBUILD_SHARED_LIBS': '1',
                  'DHYPRE_INSTALL_PREFIX': os.path.join(prefix, 'mfem'),
                  'DHYPRE_ENABLE_SHARED': '1',
                  'DCMAKE_INSTALL_PREFIX':os.path.join(prefix, 'mfem'),
                  'DCMAKE_INSTALL_NAME_DIR':os.path.join(prefix, 'mfem', 'lib'),
                  'DCMAKE_C_COMPILER': mpicc_command} 

    cmake('..', **cmake_opts)

    make('hypre')
    make_install('hypre')

    os.chdir(pwd)

def make_metis(use_int64=False, use_real64=False):
    '''
    build metis 
    '''
    path = os.path.join(extdir, 'metis')
    if not os.path.exists(path):
        assert False, "metis is not downloaded"

    pwd = chdir(path)
    
    sed_command = find_command('sed')
    if sed_command is None:
        assert False, "sed is not foudn"

    if use_int64:
        command = [sed_command, '-i',
                   's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/g',
                   'include/metis.h']
    else:
        command = [sed_command, '-i',
                   's/#define IDXTYPEWIDTH 64/#define IDXTYPEWIDTH 32/g',
                   'include/metis.h']

    if use_real64:
        command = [sed_command, '-i',
                   's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/g',
                   'include/metis.h']
    else:
        command = [sed_command, '-i',
                   's/#define REALTYPEWIDTH 64/#define REALTYPEWIDTH 32/g',
                   'include/metis.h']
    make_call(command)

    command = ['make', 'config', 'shared=1',
               'prefix='+os.path.join(prefix, 'mfem'),
               'cc='+cc_command]
    make_call(command)
    make('metis')
    make_install('metis')

    if platform == "darwin":
        command = ['install_name_tool',
                   '-id',
                   os.path.join(prefix, 'lib', 'libmetis.dylib'),
                   os.path.join(prefix, 'lib', 'libmetis.dylib'),]
        make_call(command)
    os.chdir(pwd)

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

    print("prefix here", prefix)
    
    cmake_opts = {'DCMAKE_VERBOSE_MAKEFILE':'1',
                  'DBUILD_SHARED_LIBS': '1',
                  'DMFEM_ENABLE_EXAMPLES': '1',
                  'DMFEM_ENABLE_MINIAPPS': '1',
                  'DCMAKE_SHARED_LINKER_FLAGS':'',
                  'DCMAKE_CXX_FLAGS': cxx11_flag}

    if serial:
        cmake_opts['DCMAKE_CXX_COMPILER'] = cxx_command
        cmake_opts['DMFEM_USE_EXCEPTIONS'] = '1'
        cmake_opts['DCMAKE_INSTALL_PREFIX'] = os.path.join(prefix, 'mfem', 'ser')
    else:
        cmake_opts['DCMAKE_CXX_COMPILER'] = mpicxx_command
        cmake_opts['DMFEM_USE_EXCEPTIONS'] = '0'
        cmake_opts['DCMAKE_INSTALL_PREFIX'] = os.path.join(prefix, 'mfem', 'par')
        cmake_opts['DMFEM_USE_MPI'] = '1'
        cmake_opts['DMFEM_USE_METIS_5'] = '1'
        cmake_opts['DHYPRE_DIR'] = os.path.join(hypre_prefix)
        #cmake_opts['DHYPRE_INCLUDE'] = os.path.join(hypre_prefix, 'include')
        cmake_opts['DMETIS_DIR'] = os.path.join(metis_prefix)
        #cmake_opts['DMETIS_INCLUDE'] = os.path.join(metis_prefix, 'include')
        #cmake_opts['DCMAKE_CXX_STANDARD_LIBRARIES'] = "-lHYPRE -lmetis"

        lflags = "-L" + os.path.join(metis_prefix, 'lib')
        lflags = lflags + " -L" + os.path.join(hypre_prefix, 'lib')
        cmake_opts['DCMAKE_SHARED_LINKER_FLAGS'] = lflags
        #cmake_opts['DCMAKE_EXT_LINKER_FLAGS'] = lflags

        if enable_pumi:
            cmake_opts['DMFEM_USE_PUMI'] = '1'
            cmake_opts['DPUMI_DIR'] = pumi_prefix

    pwd = chdir(path)
    cmake('..', **cmake_opts)

    txt = 'serial' if serial else 'parallel'
    make('mfem_'+txt)
    make_install('mfem_'+txt)
    
    os.chdir(pwd)

def write_setup_local():
    '''
    create setup_local.py. parameters written here will be read
    by setup.py in mfem._ser and mfem._par
    '''
    import numpy

    if build_mfem:
        mfemser = os.path.join(prefix, 'mfem', 'ser')
        mfempar = os.path.join(prefix, 'mfem', 'par')
    else:
        mfemser = mfem_prefix
        mfempar = mfem_prefix

    params = {'cxx_ser': cxx_command,
              'cc_ser': cc_command,
              'cxx_par': mpicxx_command,
              'cc_par': mpicc_command,
              'whole_archive': '--whole-archive',
              'no_whole_archive': '--no-whole-archive',
              'nocompactunwind': '',
              'swigflag': '-Wall -c++ -python -fastproxy -olddefs -keyword',

              'hypreinc': os.path.join('hypre_prefix', 'include'),
              'hyprelib': os.path.join('hypre_prefix', 'lib'),
              'metisinc': os.path.join('metis_prefix', 'include'),
              'metislib': os.path.join('metis_prefix', 'lib'),
              'numpync': numpy.get_include(),
              'mfembuilddir': os.path.join(mfempar, 'include'),
              'mfemincdir': os.path.join(mfempar, 'include', 'mfem'),
              'mfemlnkdir': os.path.join(mfempar, 'lib'),
              'mfemserbuilddir': os.path.join(mfemser, 'include'),
              'mfemserincdir': os.path.join(mfemser, 'include', 'mfem'),
              'mfemserlnkdir': os.path.join(mfemser, 'lib'),
              'add_pumi': '',
              'add_strumpack': '',
              }


    try:
        import mpi4py ## avaialbility of this is checked before
        params['mpi4pyinc'] = mpi4py.get_include()
    except ImportError:
        params['mpi4pyinc'] = ''

    def add_extra(xxx):
        params['add_' + xxx] = '1'
        params[xxx + 'inc'] = os.path.join(xxx + '_prefix', 'include')
        params[xxx + 'lib'] = os.path.join(xxx + '_prefix', 'lib')

    if enable_pumi:
        add_extra('pumi')
    if enable_strumpack:
        add_extra('strumpack')


    pwd = chdir(rootdir)
    
    fid = open('setup_local.py', 'w')
    fid.write("#  setup_local.py \n")
    fid.write("#  generated from write_setup_local.py\n")
    fid.write("#  do not edit this directly\n")
    fid.write("#  instead use Make setup_local.py\n")

    for key, value in params.items():
        text = key.lower() + ' = "' + value
        fid.write(text+"\n")
    fid.close()
    
    os.chdir(pwd)

def generate_wrapper():
    '''
    run swig.
    '''
    def ifiles():
        files = os.listdir()
        ifiles = [x for x in files if x.endswith('.i')]
        return ifiles

    def check_new(ifile):
        wfile = ifile[:-2]+'_wrap.cxx'
        if not os.path.exists(wfile):
            return True
        return os.path.getmtime(ifile) > os.path.getmtime(wfile)

    if build_mfem:
        mfemser = os.path.join(prefix, 'mfem', 'ser')
        mfempar = os.path.join(prefix, 'mfem', 'par')
    else:
        mfemser = mfem_prefix
        mfempar = mfem_prefix

    swig_command = find_command('swig')
    if swig_command is None:
        assert False, "SWIG is not installed"

    swigflag = '-Wall -c++ -python -fastproxy -olddefs -keyword'.split(' ')

    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))

    serflag = ['-I'+ os.path.join(mfemser, 'include'),
               '-I'+ os.path.join(mfemser, 'include', 'mfem')]
    for file in ifiles():
        if not check_new(file):
            continue
        command = [swig_command] + swigflag + serflag + [file]
        make_call(command)

    if not build_parallel:
        os.chdir(pwd)
        return

    chdir(os.path.join(rootdir, 'mfem', '_par'))

    import mpi4py
    parflag = ['-I'+ os.path.join(mfempar, 'include'),
               '-I'+ os.path.join(mfempar, 'include', 'mfem'),
               '-I'+ mpi4py.get_include()]

    if enable_pumi:
        parflag.append('-I'+os.path.join(pumi_prefix, 'include'))
    if enable_strumpack:
        parflag.append('-I'+os.path.join(strumpack_prefix, 'include'))

    for file in ifiles():
        if not check_new(file):
            continue
        command = [swig_command] + swigflag + parflag + [file]
        make_call(command)

    os.chdir(pwd)

def clean_wrapper():

    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))
    
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.cxx')]
    remove_files(wfiles)

    chdir(os.path.join(rootdir, 'mfem', '_par'))
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.cxx')]
    remove_files(wfiles)

    chdir(pwd)

def make_mfem_wrapper(serial=True):
    '''
    compile PyMFEM wrapper code
    '''
    write_setup_local()

    if serial:
        os.chdir(os.path.join('rootdir', 'mfem', '_ser'))
    else:
        os.chdir(os.path.join('rootdir', 'mfem', '_ser'))

    python = sys.executable
    command = [python, 'setup.py', 'build_ext', '--inplace']
    make_call(command)

class Install(_install):
    '''
    called when pyton setup.py install
    '''
    user_options = _install.user_options + [
        ('with-parallel', None, 'Installed both serial and parallel version'),
        ('mfem-prefix=', None, 'Specify locaiton of mfem' +
         'libmfem.so must exits under <mfem-prefix>/lib'),
        ('hypre-prefix=', None, 'Specify locaiton of hypre' +
         'libHYPRE.so must exits under <hypre-prefix>/lib'),
        ('metis-prefix=', None, 'Specify locaiton of metis'+
         'libmetis.so must exits under <metis-prefix>/lib'),
        ('swig', None, 'Run Swig'),
        ('ext-only', None, 'Build metis, hypre, mfem(C++) only'),
        ('skip-ext', None, 'Skip building metis, hypre, mfem(C++) only'),        

        ('CC', None, 'c compiler'),
        ('CXX', None, 'c++ compiler'),
        ('MPICC', None, 'mpic compiler'),
        ('MPICXX', None, 'mpic++ compiler'),
        ('with-metis64', None, 'use 64bit int in metis'),
        ('with-pumi', None, 'enable pumi (parallel only)'),
        ('pumi-prefix=', None, 'Specify locaiton of pumi'),
        ('with-strumpack', None, 'enable strumpack (parallel only)'),
        ('strumpack-prefix=', None, 'Specify locaiton of strumpack'),
    ]

    def initialize_options(self):
        _install.initialize_options(self)

        self.swig = False
        self.ext_only = False
        self.with_parallel = False
        self.mfem_prefix = ''
        self.metis_prefix = ''
        self.hypre_prefix = ''

        self.with_metis64 = False
        self.with_pumi = False
        self.pumi_prefix = ''

        self.with_strumpack = False
        self.strumpack_prefix = ''

        self.CC = ''
        self.CXX = ''
        self.MPICC = ''
        self.MPICXX = ''        

    def run(self):
        global prefix, dry_run, verbose
        global build_mfem, build_parallel, build_serial
        global mfem_prefix, metis_prefix, hypre_prefix
        global cc_command, cxx_command, mpicc_command, mpicxx_command
        global enable_pumi, pumi_prefix
        global enable_strumpack, strumpack_prefix

        prefix = os.path.expanduser(self.prefix)
        dry_run = self.dry_run
        verbose = bool(self.verbose)
        
        swig_only = bool(self.swig)
        ext_only = bool(self.ext_only)
        build_wrapper = (not swig_only and not ext_only)

        mfem_prefix = os.path.expanduser(self.mfem_prefix)
        build_parallel = bool(self.with_parallel)     # controlls PyMFEM parallel
        metis_64 = bool(self.with_metis64)
        enable_pumi = bool(self.with_pumi)
        enable_strumpack = bool(self.with_strumpack)


        if build_parallel:
            try:
                import mpi4py
            except ImportError:
                assert False, "Can not import mpi4py"

        if mfem_prefix != '':
            path = os.path.join(mfem_prefix, 'lib', 'libmfem'+dylibext)
            assert os.path.exists(path), "libmfem.so is not found in the specified <path>/lib"
            build_mfem = False
            hypre_prefix = mfem_prefix
            metis_prefix = mfem_prefix
        else:
            build_mfem = True
            build_hypre = build_parallel
            build_metis = build_parallel
            hypre_prefix = os.path.join(prefix, 'mfem')
            metis_prefix = os.path.join(prefix, 'mfem')

        if self.hypre_prefix != '':
            hypre_prefix = os.path.expanduser(self.hypre_prefix)            
            path = os.path.join(hypre_prefix, 'lib', 'libHYPRE'+dylibext)
            assert os.path.exists(path), "libHYPRE.so is not found in the specified <path>/lib"
            build_hypre = False

        if self.metis_prefix != '':
            metis_prefix = os.path.expanduser(self.metis_prefix)
            path = os.path.join(metis_prefix, 'lib', 'libmetis'+dylibext)
            assert os.path.exists(path), "libmetis.so is not found in the specified <path>/lib"
            build_metis = False

        if self.pumi_prefix != '':
            pumi_prefix = os.path.expanduser(self.pumi_prefix)
        else:
            pumi_prefix = mfem_prefix

        if self.strumpack_prefix != '':
            strumpack_prefix = os.path.expanduser(self.strumpack_prefix)
        else:
            strumpack_prefix = mfem_prefix

        if self.CC != '':
            cc_command = self.CC
        if self.CXX != '':
            cxx_command = self.CXX
        if self.MPICC != '':
            mpicc_command = self.MPICC
        if self.MPICXX != '':
            mpicxx_command = self.MPICXX

            
        print("----configuration----")
        print(" prefix", prefix)
        print(" when needed, the dependency (mfem/hypre/metis) will be installed under " +
              prefix + "/mfem")
        print(" build mfem : " + ("Yes" if build_mfem else "No"))        
        print(" build metis : " + ("Yes" if build_metis else "No"))
        print(" build hypre : " + ("Yes" if build_hypre else "No"))
        print(" build serial wrapper: Yes")
        print(" build parallel wrapper : " + ("Yes" if build_parallel else "No"))
        print(" c compiler : " + cc_command)
        print(" c++ compiler : " + cxx_command)
        print(" mpi-c compiler : " + mpicc_command)
        print(" mpi-c++ compiler : " + mpicxx_command)        
        print("")

        if swig_only:
            generate_wrapper()

        else:
            
            if build_metis:
                download('metis')
                make_metis(use_int64=metis_64)
            if build_hypre:
                download('hypre')
                cmake_make_hypre()         
            if build_mfem:
                download('mfem')
                build_serial = True

                cmake_make_mfem(serial=True)
                if build_parallel:
                    cmake_make_mfem(serial=False)
                    
        if build_wrapper:                    
            _install.run(self)

class BuildPy(_build_py):
    '''
    Called when python setup.py build_py
    '''
    user_options = _build_py.user_options
    def initialize_options(self):
        _build_py.initialize_options(self)
        self.custom_option = None
        
    def finalize_options(self):
        _build_py.finalize_options(self)
        
    def run(self):
        if not build_mfem:
            clean_wrapper()
            generate_wrapper()

        make_mfem_wrapper(serial=True)
        if build_parallel:
            make_mfem_wrapper(serial=False)
        
        _build_py.run(self)
    
class Clean(_clean):
    '''
    Called when python setup.py clean
    '''
    user_options = _clean.user_options + [
        ('externals', None, 'clean exteranal dependencies)'),
        ('mfem', None, 'clean mfem'),
        ('metis', None, 'clean metis'),
        ('hypre', None, 'clean hypre'),
        ('swig', None,  'clean swig'),        
        ('all-externals', None, 'delete all externals'),
        ]

    def initialize_options(self):
        _clean.initialize_options(self)
        self.externals = False
        self.mfem = False
        self.hypre = False
        self.metis = False
        self.swig = False
        self.all_externals = False

    def run(self):
        global dry_run, verbose
        dry_run = self.dry_run
        verbose = bool(self.verbose)
        
        os.chdir(extdir)

        make_command = find_command('make')
        
        if self.externals or self.mfem:
            path = os.path.join(extdir, 'mfem', 'cmbuild_par')
            if os.path.exists(path):
                shutil.rmtree(path)
            path = os.path.join(extdir, 'mfem', 'cmbuild_ser')
            if os.path.exists(path):
                shutil.rmtree(path)
        if self.externals or self.hypre:    
            path = os.path.join(extdir, 'hypre', 'cmbuild')
            if os.path.exists(path):
                shutil.rmtree(path)
        if self.externals or self.metis:    
            path = os.path.join(extdir, 'metis')
            if os.path.exits(path):
                os,chdir(path)
                command = ['make', 'clean']
                subprocess.check_call(command)                
        if self.all_externals or self.hypre:
            for xxx in ('metis', 'hypre', 'mfem'):
                path = os.path.join(extdir, xxx)
                if os.path.exists(path):
                    shutil.rmtree(path)
        if self.swig:
            clean_wrapper()
                      
        os.chdir(rootdir)
        _clean.run(self)

datafiles = [os.path.join('data', f) for f in os.listdir('data')]
def run_setup():
    setup_args = metadata.copy()
    
    setup(
        cmdclass = {'build_py': BuildPy,
                    'install': Install,
                    'clean': Clean},
        install_requires=["numpy"],
        packages=find_packages(),
        extras_require={},
        package_data={'mfem._par':['*.so'], 'mfem._ser':['*.so']},
        data_files=[('data', datafiles)],
        entry_points={},
        **setup_args)

def main():
    run_setup()

if __name__ == '__main__':
    main()
