#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
print('building paralel version')

## first load variables from PyMFEM_ROOT/setup_local.py
import sys
import os


ddd = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
root =  os.path.abspath(os.path.join(ddd, '..', '..'))

sys.path.insert(0, root)
from  setup_local import *

## this forces to use compiler written in setup_local.py
if cc_par != '': os.environ['CC'] = cc_par
if cxx_par != '': os.environ['CXX'] = cxx_par

from distutils.core import setup, Extension
from distutils.core import *
from distutils      import sysconfig

modules= ["io_stream", "vtk", "sort_pairs", "datacollection",
          "globals", "mem_manager", "device", "hash", "stable3d",
          "cpointers",
          "error", "array", "common_functions",
          "segment", "point", "hexahedron", "quadrilateral",
          "tetrahedron", "triangle", "wedge",
          "socketstream", "handle",
          "blockvector", "blockoperator", "blockmatrix",
          "vertex", "sets", "element", "table",
          "fe", "mesh", "fespace",
          "fe_coll", "coefficient",
          "linearform", "vector", "lininteg", "complex_operator",
          "complex_fem",          
          "gridfunc", "hybridization", "bilinearform",
          "bilininteg", "intrules", "sparsemat", "densemat",
          "solvers", "estimators", "mesh_operators", "ode",
          "sparsesmoothers", "ncmesh",
          "matrix", "operators", "eltrans", "geom",
          "nonlininteg", "nonlinearform",
          "pmesh", "pncmesh", "communication",
          "pfespace", "pgridfunc",
          "plinearform", "pbilinearform", "pnonlinearform",
          "hypre", "restriction", "prestriction",
          "fespacehierarchy", "multigrid", "constraints"]

if add_pumi != '':
    modules.append("pumi")
extra_compile_args = [cxx11flag, '-DSWIG_TYPE_TABLE=PyMFEM']

sources = {name: [name + "_wrap.cxx"] for name in modules}

proxy_names = {name: '_'+name for name in modules}

import numpy
numpyinc = numpy.get_include()
import mpi4py
mpi4pyinc = mpi4py.get_include()

libraries    = ['mfem', 'HYPRE', 'metis']
include_dirs = [mfembuilddir, mfemincdir, numpyinc, mpi4pyinc, hypreinc, metisinc]
#                mpichinc, hypreinc,]
library_dirs = [mfemlnkdir, hyprelib, metis5lib,]

if add_pumi != '':
   include_dirs.append(pumiinc)
   library_dirs.append(pumilib)
   
if add_strumpack:
    modules.append("strumpack")
    extra_compile_args.append('-std=c++11')
    sources["strumpack"] = ["strumpack_wrap.cxx"]
    proxy_names["strumpack"] = "_strumpack"
    if strumpackinc != "":
        include_dirs.append(strumpackinc)
    if strumpacklib != "":
        library_dirs.append(strumpacklib)
        
if add_cuda:
    include_dirs.append(cudainc)

if add_libceed:
    include_dirs.append(libceedinc)

import six
if six.PY3:
    macros = [('TARGET_PY3', '1'),]
else:
    macros = []
        
ext_modules = [Extension(proxy_names[modules[0]],
                        sources=sources[modules[0]],
                        extra_compile_args = extra_compile_args,
                        extra_link_args = [],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        runtime_library_dirs = library_dirs,
                         libraries = libraries,
                         define_macros=macros),]

ext_modules.extend([Extension(proxy_names[name],
                        sources=sources[name],
                        extra_compile_args = extra_compile_args,
                        extra_link_args = [],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        runtime_library_dirs = library_dirs,
                        libraries = libraries,
                        define_macros=macros)
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

setup (name = 'mfem_parallel',
       version = version,
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
