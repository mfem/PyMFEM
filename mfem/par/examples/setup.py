#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
from distutils.core import setup, Extension
import os
from os.path import expanduser

from distutils.core import *
from distutils      import sysconfig

import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
     numpy_include = numpy.get_numpy_include()

home = expanduser("~")
mfem_dir = os.path.join(home, 'mfem-3.1')+'/'
mfem_lib ='mfem'

ex3_module = Extension('_ex3',
                        sources=['ex3_wrap.cxx'],
                        include_dirs=[mfem_dir],
                        library_dirs=[mfem_dir],
                        libraries=[mfem_lib])
waveguide_module = Extension('_waveguide',
                        sources=['waveguide_wrap.cxx'],
                        include_dirs=[mfem_dir],
                        library_dirs=[mfem_dir],
                        libraries=[mfem_lib])

setup (name = 'mfem',
       version = '0.1',
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = [ex3_module, waveguide_module],
       py_modules = ["ex3", "waveguide"], 
       )
