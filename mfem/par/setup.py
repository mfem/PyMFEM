#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
print('building paralel version')

## first load variables from PyMFEM_ROOT/setup_local.py
import sys
import os
root =  os.path.abspath(os.path.join(os.path.realpath(__file__),
                                     '..', '..', '..'))
sys.path.insert(0, root)
from  setup_local import *

## this forces to use compiler written in setup_local.py
os.environ['CXX'] = mpicxx
os.environ['CC'] = mpicc

from distutils.core import setup, Extension
from distutils.core import *
from distutils      import sysconfig

modules= [
          "array",
          "blockvector", "blockoperator", "blockmatrix",
          "vertex", "sets", "element", "table",
          "fe", "mesh", "fespace", 
          "fe_coll", "coefficient",
          "linearform", "vector", "lininteg",
          "gridfunc", "bilinearform",
          "bilininteg", "intrules", "sparsemat", "densemat",
          "solvers", "estimators", "mesh_operators", "ode",
          "sparsesmoothers",
          "matrix", "operators", "ncmesh", "eltrans", "geom",
          "nonlininteg", "nonlinearform", 
          "pmesh", "pncmesh", "communication",
          "pfespace", "pgridfunc",
          "plinearform", "pbilinearform", "pnonlinearform",
          "hypre"]

sources = {name: [name + "_wrap.cxx"] for name in modules}
#sources['solvers'] = ['solvers_p_wrap.cxx']

proxy_names = {name: '_'+name for name in modules}

extra_link_args0 =  ['-Wl,'+whole_archive+','+mfemlnkdir+'/libmfem.a'+no_whole_archive]
extra_link_args =  [metis4liba]
include_dirs = [mfemincdir, numpyincdir, mpi4pyincdir,
                mpichincdir, hypreincdir]
library_dirs = [mfemlnkdir, hyprelnkdir]
libraries    = [hyprelib]

ext_modules = [Extension(proxy_names[modules[0]],
                        sources=sources[modules[0]],
                        extra_link_args = extra_link_args + extra_link_args0,
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                         libraries = libraries)]

ext_modules.extend([Extension(proxy_names[name],
                        sources=sources[name],
                        extra_link_args = extra_link_args,
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                              libraries = libraries)
               for name in modules[1:]])

#ext_modules = [Extension(proxy_names[name],
#                        sources=sources[name],
#                        extra_link_args = extra_link_args,
#                        include_dirs = include_dirs,
#                        library_dirs = library_dirs,
#                        libraries = libraries)
#              for name in modules]


setup (name = 'mfem',
       version = '0.2',
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
