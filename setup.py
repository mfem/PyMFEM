"""

  MFEM setuptools based setup module.
  

  python setup.py build  # build mfem and PyMFEM in serial
  python setup.py build --parallel # build metis/hypre/mfem and PyMFEM in parallel


  (plan)
  python setup.py build --cuda


  (note we will install evertying under <prefix>/mfem/)
     <prefix> /mfem/par  : mfem parallel
     <prefix> /mfem/ser  : mfem serial 

     <package location>/sandbox/TwoPi : package manager
     <package location>/sandbox/src : download area

"""
import os
from os import path

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.install import install as _install

import subprocess
# To use a consistent encoding
from codecs import open

import urllib
import gzip
import re


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

repos = {"mfem": "https://github.com/mfem/mfem/archive/v4.2.tar.gz",
         "metis": "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz",
         "hypre": "https://github.com/hypre-space/hypre/archive/v2.18.2.tar.gz",}

rootdir = path.abspath(path.dirname(__file__))
extdir = path.join(rootdir, 'external')

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

def download(url):
    from urllib import request
    import tarfile
    ftpstream = request.urlopen(url)
    targz = tarfile.open(fileobj=ftpstream, mode="r|gz")
    targz.extractall()
    return targz.getnames()[0].split('/')[0]

class Install(_install):
    user_options = _install.user_options + [
        ('with-parallel', None, 'Installed both serial and parallel version'),
        ('mfem-prefix=', None, 'MFEM prefix if when is already is installed.'+
                               'Needs to be dynamic library with specific version'),
    ]

    def initialize_options(self):
        _install.initialize_options(self)
        self.with_metis = False
        self.with_hypre = False
        self.with_parallel = False
        self.mfem_prefix = ''
        #self.someval = None

    def finalize_options(self):
        #print("vlue of someopt is", self.someopt)
        _install.finalize_options(self)

    def run(self):
        global with_metis, with_hypre, with_parallel, mfem_prefix

        mfem_prefix = self.mfem_prefix
        if mfem_prefix != '':
            path = os.path.join(mfem_prefix, 'lib', 'libmfem.so')
            assert os.path.exists(path), "libmfem.so is not found in the specified <path>/lib"

        with_parallel = bool(self.with_parallel)
        with_hypre = with_parallel and mfem_prefix != ''
        with_metis = with_parallel and mfem_prefix != ''

        print("configuration")
        print(" build serial : Yes")
        print(" build parallel : " + ("Yes" if with_parallel else "No"))
        print(" build metis : " + ("Yes" if with_metis else "No"))
        print(" build hypre : " + ("Yes" if with_hypre else "No"))
        print(" MFEM to be linked : " + ("default" if mfem_prefix == '' else "under "+mfem_prefix))
        
class BuildPy(_build_py):
    user_options = _build_py.user_options + [
                   ('custom-option=', None, 'Path to something')
                   ]

    def initialize_options(self):
        print(build_py.user_options)
        _build_py.initialize_options(self)
        self.custom_option = None
        
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


datafiles = [path.join('data', f) for f in os.listdir('data')]
def run_setup():
    setup_args = metadata.copy()
    
    setup(
        cmdclass = {'build_py': BuildPy,
                    'install': Install},
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
