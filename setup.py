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
import os
from os import path

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.install import install as _install
try:
    from setuptools._distutils.command.clean import clean as _clean
except ImportError:
    from distutils.command.clean import clean as _clean

import subprocess
from subprocess import DEVNULL 
# To use a consistent encoding
from codecs import open

import urllib
import gzip
import re
import shutil

### constants
repos = {"mfem": "https://github.com/mfem/mfem/archive/v4.2.tar.gz",
         "metis": "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz",
         "hypre": "https://github.com/hypre-space/hypre/archive/v2.18.2.tar.gz",}

rootdir = path.abspath(path.dirname(__file__))
extdir = path.join(rootdir, 'external')

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
    VERSIONFILE = path.join('mfem', '__init__.py')
    initfile_lines = open(VERSIONFILE, 'rt').readlines()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    for line in initfile_lines:
        mo = re.search(VSRE, line, re.M)
        if mo:
            return mo.group(1)
    raise RuntimeError('Unable to find version string in %s.' % (VERSIONFILE,))

def long_description():
    with open(path.join(rootdir, 'README')) as f:
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

def make_call(command):
    '''
    call command
    '''
    kwargs = {}
    if not verbose:
        kwargs['stdout'] = DEVNULL
        kwargs['stderr'] = DEVNULL
    subprocess.check_call(command, **kwargs)

def make(target):
    import multiprocessing
    command = ['make', '-j', str(max((multiprocessing.cpu_count()-1, 1)))]

    try:
        if dry_run:
            print("calling ... " + " ".join(command))
        else:
            make_call(command)
    except subprocess.CalledProcessError:
        assert False, "Error in building: " + target

def make_install(target):
    command = ['make', 'install']
    try:
        if dry_run:
            print("calling ... " + " ".join(command))
        else:
            make_call(command)
    except subprocess.CalledProcessError:
        assert False, "Failed when calling install: " + target

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
    
    if dry_run:
        return
        
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

    try:
        if dry_run:
            print("calling ... " + " ".join(command))
        else:
            make_call(command)            
    except:
        assert False, "Failed to call cmake"

def cmake_make_hypre():
    pass

def make_hypre():
    pass

def cmake_make_mfem(serial=True):
    cmbuild = 'cmbuild_ser' if serial else 'cmbuild_par'
    path = os.path.join(extdir, 'mfem', cmbuild)
    if os.path.exists(path):
        print("working directory already exists!") 
    else:
        os.makedirs(path)
    
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
        cmake_opts['DHYPRE_DIR'] = os.path.join(hypre_prefix, 'lib')
        cmake_opts['DHYPRE_INCLUDE'] = os.path.join(hypre_prefix, 'include')
        cmake_opts['DMETIS_DIR'] = os.path.join(metis_prefix, 'lib')
        cmake_opts['DMETIS_INCLUDE'] = os.path.join(metis_prefix, 'include')
        cmake_opts['DCMAKE_CXX_STANDARD_LIBRARIES']="-lHYPRE -lmetis"
        cmake_opts['DCMAKE_SHARED_LINKER_FLAGS'] = "-L" + metis_prefix + " -L" + hypre_prefix
        cmake_opts['DCMAKE_EXT_LINKER_FLAGS'] = "-L"+metis_prefix + " -L"+hypre_prefix
        
        if enable_pumi:
            cmake_opts['DMFEM_USE_PUMI'] = '1'
            cmake_opts['DPUMI_DIR'] = pumi_prefix

    os.chdir(path)
    cmake('..', **cmake_opts)

    txt = 'serial' if serial else 'parallel'
    make('mfem_'+txt)
    make_install('mfem_'+txt)
    
class Install(_install):
    user_options = _install.user_options + [
        ('parallel', None, 'Installed both serial and parallel version'),
        ('mfem-prefix=', None, 'Specify locaiton of mfem' + 
                                'libmfem.so must exits under <mfem-prefix>/lib'),
        ('hypre-prefix=', None, 'Specify locaiton of hypre' + 
                                'libHYPRE.so must exits under <hypre-prefix>/lib'),
        ('metis-prefix=', None, 'Specify locaiton of metis'+
                                'libmetis.so must exits under <metis-prefix>/lib'),
        ('pumi-prefix=', None, 'Specify locaiton of pumi'),
        ('CC', None, 'c compiler'),
        ('CXX', None, 'c++ compiler'),
        ('MPICC', None, 'mpic compiler'),
        ('MPICXX', None, 'mpic++ compiler'),
        ('with-pumi', None, 'enable pumi (parallel only)'),
    ]

    def initialize_options(self):
        _install.initialize_options(self)
        self.parallel = False
        self.mfem_prefix = ''
        self.metis_prefix = ''
        self.hypre_prefix = ''
        
        self.with_pumi = False
        self.pumi_prefix = ''
        
        self.CC = ''
        self.CXX = ''
        self.MPICC = ''
        self.MPICXX = ''        
        #self.someval = None

    def run(self):
        global prefix, dry_run, verbose
        global build_parallel, build_serial
        global mfem_prefix, metis_prefix, hypre_prefix
        global cc_command, cxx_command, mpicc_command, mpicxx_command
        global enable_pumi, pumi_prefix

        prefix = self.prefix
        dry_run = self.dry_run
        
        mfem_prefix = self.mfem_prefix
        build_parallel = bool(self.parallel)     # controlls PyMFEM parallel
        enable_pumi = bool(self.with_pumi)
        verbose = bool(self.verbose)

        if mfem_prefix != '':
            path = os.path.join(mfem_prefix, 'lib', 'libmfem'+dylibext)
            assert os.path.exists(path), "libmfem.so is not found in the specified <path>/lib"
            build_mfem = False
            hypre_prefix = self.mfem_prefix
            metis_prefix = self.mfem_prefix
        else:
            build_mfem = True            
            build_hypre = build_parallel
            build_metis = build_parallel
            hypre_prefix = self.prefix
            metis_prefix = self.prefix

        if self.hypre_prefix != '':
            path = os.path.join(self.hypre_prefix, 'lib', 'libHYPRE'+dylibext)
            assert os.path.exists(path), "libHYPRE.so is not found in the specified <path>/lib"
            build_hypre = False
            hypre_prefix = self.hypre_prefix

        if self.metis_prefix != '':
            path = os.path.join(self.metis_prefix, 'lib', 'libmetis'+dylibext)
            assert os.path.exists(path), "libmetis.so is not found in the specified <path>/lib"
            build_metis = False
            metis_prefix = self.metis_prefix

        if self.pumi_prefix != '':
            pumi_prefix = self.pumi_prefix
        else:
            pumi_prefix = self.mfem_prefix

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
        
        if build_metis:
            download('metis')
            make_metis()
        if build_hypre:
            download('hypre')
            cmake_make_hypre()            
        if build_mfem:
            download('mfem')
            cmake_make_mfem(serial=True)
            if build_parallel:
                cmake_make_mfem(serial=False)

        _install.run(self)
                
class BuildPy(_build_py):
    user_options = _build_py.user_options + [
                   ('custom-option=', None, 'Path to something')
                   ]

    def initialize_options(self):
        _build_py.initialize_options(self)
        self.custom_option = None
        
    def finalize_options(self):
        _build_py.finalize_options(self)
        
    def run(self):
        assert False, "stop here"
        
        os.chdir(extdir)
        for repo in repos:
            try:
                print("downloading", repo, " from ", repos[repo] )
                download(repos[repo])
            except BaseException:
                assert False, "Failed to download dependency:" + repo
                
        os.chdir(rootdir)
        _build_py.run(self)

class Clean(_clean):
    user_options = _clean.user_options + [
                   ('externals', None, 'clean exteranal dependencies)'),
                   ('mfem', None, 'clean mfem'),
                   ('metis', None, 'clean metis'),
                   ('hypre', None, 'clean hypre'),
                   ('all-externals', None, 'delete all externals'),
                   ]

    def initialize_options(self):
        _clean.initialize_options(self)
        self.externals = False
        self.mfem = False
        self.hypre = False
        self.metis = False
        self.all_externals = False

    def run(self):
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
        os.chdir(rootdir)
        _clean.run(self)

datafiles = [path.join('data', f) for f in os.listdir('data')]
def run_setup():
    setup_args = metadata.copy()
    
    setup(
        cmdclass = {'build_py': BuildPy,
                    'install': Install,
                    'clean': Clean},
        install_requires=[],
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
