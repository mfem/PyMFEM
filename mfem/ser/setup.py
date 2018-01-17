#!/usr/bin/env python

"""
Serial version setup file
"""
print('building serial version')

## first load variables from PyMFEM_ROOT/setup_local.py
import sys
import os
root =  os.path.abspath(os.path.join(os.path.realpath(__file__),
                                     '..', '..', '..'))
sys.path.insert(0, root)
from  setup_local import *

## this forces to use compiler written in setup_local.py
if cc_ser != '': os.environ['CC'] = cc_ser
if cxx_ser != '': os.environ['CXX'] = cxx_ser

from distutils.core import *
from distutils      import sysconfig

modules= ["array", "common_functions", "socketstream",
          "blockvector", "blockoperator", "blockmatrix",
          "vertex", "sets", "element", "table", "fe",
          "mesh", "fespace", 
          "fe_coll", "coefficient",
          "linearform", "vector", "lininteg",
          "gridfunc", "hybridization", "bilinearform",
          "bilininteg", "intrules", "sparsemat", "densemat",
          "solvers", "estimators", "mesh_operators", "ode",
          "sparsesmoothers",
          "matrix", "operators", "ncmesh", "eltrans", "geom",
          "nonlininteg", "nonlinearform",
          "ostream_typemap", "istream_typemap"]

sources = {name: [name + "_wrap.cxx"] for name in modules}
proxy_names = {name: '_'+name for name in modules}

#extra_text = [x for x in
#              ['-Wl', whole_archive,  mfemserlnkdir+'/libmfem.a',
#               no_whole_archive] if x != '']
#extra_link_args =  [','.join(extra_text)]
extra_link_args =  []
include_dirs = [mfemserbuilddir, mfemserincdir, numpyincdir, boostincdir]
library_dirs = [mfemserlnkdir, boostlibdir]
libraries    = [mfemlib, boostlib]


ext_modules = [Extension(proxy_names[modules[0]],
                         sources=sources[modules[0]],
                         extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],
                         extra_link_args = extra_link_args,
                         include_dirs = include_dirs,
                         library_dirs = library_dirs,
                         runtime_library_dirs = library_dirs,
                         libraries = libraries)]

ext_modules.extend([Extension(proxy_names[name],
                              sources=sources[name],
                              extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],                              
                              extra_link_args = [],
                              include_dirs = include_dirs,
                              library_dirs = library_dirs,
                              libraries = libraries)
               for name in modules[1:]])


setup (name = 'mfem',
       version = '0.1',
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
