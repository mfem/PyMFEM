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
if cc_par != '': os.environ['CC'] = cc_par
if cxx_par != '': os.environ['CXX'] = cxx_par

from distutils.core import setup, Extension
from distutils.core import *
from distutils      import sysconfig

modules= ["cpointers",
          "error", "array", "common_functions",
          "point", "segment", 
          "socketstream", "handle", 
          "blockvector", "blockoperator", "blockmatrix",
          "vertex", "sets", "element", "table",
          "fe", "mesh", "fespace", 
          "fe_coll", "coefficient",
          "linearform", "vector", "lininteg",
          "gridfunc", "hybridization", "bilinearform",
          "bilininteg", "intrules", "sparsemat", "densemat",
          "solvers", "estimators", "mesh_operators", "ode",
          "sparsesmoothers", "ncmesh", 
          "matrix", "operators", "eltrans", "geom",
          "nonlininteg", "nonlinearform", 
          "pmesh", "pncmesh", "communication",
          "pfespace", "pgridfunc",
          "plinearform", "pbilinearform", "pnonlinearform",
          "hypre", 
          "ostream_typemap", "istream_typemap"]

sources = {name: [name + "_wrap.cxx"] for name in modules}
#sources['solvers'] = ['solvers_p_wrap.cxx']

proxy_names = {name: '_'+name for name in modules}


#extra_text = [x for x in
#              ['-Wl', whole_archive, metisliba, mfemlnkdir+'/libmfem.a',
#               no_whole_archive] if x != '']
#extra_link_args0 = []
#extra_link_args =  []
#libraries0 = [] if metisliba != '' else [metislib]
#libraries0.extend([hyprelib, boostlib])
libraries    = [libboostiostreams, 'HYPRE', 'mfem']

include_dirs = [mfembuilddir, mfemincdir, numpyinc, mpi4pyinc,
                mpichinc, hypreinc, boostinc]
library_dirs = [mfemlnkdir, hyprelib, metis5lib, boostlib]


ext_modules = [Extension(proxy_names[modules[0]],
                        sources=sources[modules[0]],
                        extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],                            
                        extra_link_args = [],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        runtime_library_dirs = library_dirs,  
                        libraries = ['metis']+libraries)]

ext_modules.extend([Extension(proxy_names[name],
                        sources=sources[name],
                        extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'], 
                        extra_link_args = [],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        libraries = libraries)
               for name in modules[1:]])

setup (name = 'mfem_parallel',
       version = '3.3.3',
       author      = "S.Shiraiwa",
       description = """MFEM wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
