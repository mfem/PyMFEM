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

from distutils.core import *
from distutils      import sysconfig

modules= ["array",
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
          "nonlininteg", "nonlinearform", ]

sources = {name: [name + "_wrap.cxx"] for name in modules}
proxy_names = {name: '_'+name for name in modules}

extra_link_args =  ['-Wl,'+whole_archive+','+mfemserlnkdir+'/libmfem.a']
include_dirs = [mfemserincdir, numpyincdir]
library_dirs = [mfemserlnkdir]
libraries    = [mfemserlib]


ext_modules = [Extension(proxy_names[modules[0]],
                        sources=sources[modules[0]],
                        extra_link_args = extra_link_args,
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                         libraries = [])]

ext_modules.extend([Extension(proxy_names[name],
                        sources=sources[name],
                        extra_link_args = [],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        libraries = [])
               for name in modules[1:]])


setup (name = 'mfem',
       version = '0.1',
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
