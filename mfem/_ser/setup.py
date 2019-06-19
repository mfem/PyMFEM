#!/usr/bin/env python

"""
Serial version setup file
"""
print('building serial version')

## first load variables from PyMFEM_ROOT/setup_local.py
import sys
import os

ddd = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
root =  os.path.abspath(os.path.join(ddd, '..', '..'))

sys.path.insert(0, root)
from  setup_local import *

## this forces to use compiler written in setup_local.py
if cc_ser != '': os.environ['CC'] = cc_ser
if cxx_ser != '': os.environ['CXX'] = cxx_ser

from distutils.core import *
from distutils      import sysconfig

modules= ["globals", "mem_manager", "device", "hash",
          "error", "array", "common_functions", "socketstream", "handle",
          "segment", "point",
          "blockvector", "blockoperator", "blockmatrix",
          "vertex", "sets", "element", "table", "fe",
          "mesh", "fespace", 
          "fe_coll", "coefficient",
          "linearform", "vector", "lininteg", "complex_operator",
          "gridfunc", "hybridization", "bilinearform",
          "bilininteg", "intrules", "sparsemat", "densemat",
          "solvers", "estimators", "mesh_operators", "ode",
          "sparsesmoothers",
          "matrix", "operators", "ncmesh", "eltrans", "geom",
          "nonlininteg", "nonlinearform",
          "ostream_typemap", "istream_typemap"]

sources = {name: [name + "_wrap.cxx"] for name in modules}
proxy_names = {name: '_'+name for name in modules}

include_dirs = [mfemserbuilddir, mfemserincdir, numpyinc, boostinc]
library_dirs = [mfemserlnkdir, boostlib]
libraries    = [libboostiostreams, 'mfem']

extra_compile_args = [cxx11flag, '-DSWIG_TYPE_TABLE=PyMFEM']

ext_modules = [Extension(proxy_names[modules[0]],
                         sources=sources[modules[0]],
                         extra_compile_args = extra_compile_args,
                         extra_link_args = [],
                         include_dirs = include_dirs,
                         library_dirs = library_dirs,
                         runtime_library_dirs = library_dirs,
                         libraries = libraries)]

ext_modules.extend([Extension(proxy_names[name],
                              sources=sources[name],
                              extra_compile_args = extra_compile_args,                              
                              extra_link_args = [],
                              include_dirs = include_dirs,
                              runtime_library_dirs = library_dirs,                              
                              library_dirs = library_dirs,
                              libraries = libraries)
               for name in modules[1:]])

### read version number from __init__.py
path = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(__file__))),
                    '..', '__init__.py')
fid = open(path, 'r')
lines=fid.readlines()
fid.close()
for x in lines:
    if x.strip().startswith('__version__'):
        version = eval(x.split('=')[-1].strip())

setup (name = 'mfem_serial',
       version = version,
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
