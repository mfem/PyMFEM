"""

  Petra-M setuptools based setup module.

"""
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path, listdir
import re

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README')) as f:
    long_description = f.read()
    
def get_version():
    VERSIONFILE = path.join('mfem', '__init__.py')
    initfile_lines = open(VERSIONFILE, 'rt').readlines()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    for line in initfile_lines:
        mo = re.search(VSRE, line, re.M)
        if mo:
            return mo.group(1)
    raise RuntimeError('Unable to find version string in %s.' % (VERSIONFILE,))

datafiles = [path.join('data', f) for f in listdir('data')]
setup(
    name='PyMFEM',
    version=get_version(),

    description='PyMFEM',
    long_description=long_description,
    url='https://github.com/mfem/PyMFEM',
    author='S. Sihraiwa',
    author_email='shiraiwa@psfc.mit.edu',
    license='LGPL-2.1',

    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics'
        'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
        'Programming Language :: Python :: 2.7',
    ],

    keywords='MFEM physics',
    packages=find_packages(),
    install_requires=[],
    extras_require={},
    package_data={'mfem.par':['*.so'], 'mfem.ser':['*.so']},
    data_files=[('data', datafiles)],
    entry_points={},
)
