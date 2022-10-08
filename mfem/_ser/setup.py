#!/usr/bin/env python

"""
Serial version setup file
"""


## first load variables from PyMFEM_ROOT/setup_local.py
import sys
import os
import numpy

from distutils.core import Extension, setup


ddd = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
root =  os.path.abspath(os.path.join(ddd, '..', '..'))

def get_version():
    ### read version number from __init__.py    
    path = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(__file__))),
                    '..', '__init__.py')
    fid = open(path, 'r')
    lines=fid.readlines()
    fid.close()
    for x in lines:
        if x.strip().startswith('__version__'):
            version = eval(x.split('=')[-1].strip())
    return version

def get_extensions():
    sys.path.insert(0, root)
    try:
        from setup_local import (mfemserbuilddir, mfemserincdir, mfemsrcdir, mfemserlnkdir,
                                 mfemstpl,
                                 cc_ser, cxx_ser,
                                 cxx11flag,                                 
                                 add_cuda, add_libceed, add_suitesparse, add_gslibs,)

        include_dirs = [mfemserbuilddir, mfemserincdir, mfemsrcdir,
                        numpy.get_include()]
        library_dirs = [mfemserlnkdir,]
        
    except ImportError:
        if 'clean' not in sys.argv:
            raise
        
        cc_ser = ''
        cxx_ser = ''
        include_dirs = []
        library_dirs = []
        mfemstpl = ''
        add_cuda = ''
        add_libceed = ''
        add_suitesparse = ''
        add_gslibs = ''
        cxx11flag = ''
        
    libraries    = ['mfem']
    
    ## remove current directory from path
    print("__file__", os.path.abspath(__file__))
    if '' in sys.path:
        sys.path.remove('')
    items = [x for x in sys.path if os.path.abspath(x) == os.path.dirname(os.path.abspath(__file__))]
    for x in items:
        sys.path.remove(x)
    print("sys path", sys.path)

    ## this forces to use compiler written in setup_local.py
    if cc_ser != '': os.environ['CC'] = cc_ser
    if cxx_ser != '': os.environ['CXX'] = cxx_ser

    modules= ["io_stream", "vtk", "sort_pairs", "datacollection",
              "cpointers", "symmat",
              "globals", "mem_manager", "device", "hash", "stable3d",
              "error", "array", "common_functions", "socketstream", "handle",
              "fe_base", "fe_fixed_order", "fe_h1", "fe_l2",
              "fe_nd", "fe_nurbs", "fe_pos", "fe_rt", "fe_ser", "doftrans",
              "segment", "point", "hexahedron", "quadrilateral",
              "tetrahedron", "triangle", "wedge",
              "blockvector", "blockoperator", "blockmatrix",
              "vertex", "sets", "element", "table", "fe",
              "mesh", "fespace",
              "fe_coll", "coefficient",
              "linearform", "vector", "lininteg", "complex_operator",
              "complex_fem",
              "gridfunc", "hybridization", "bilinearform",
              "bilininteg", "intrules", "sparsemat", "densemat",
              "solvers", "estimators", "mesh_operators", "ode",
              "sparsesmoothers",
              "matrix", "operators", "ncmesh", "eltrans", "geom",
              "nonlininteg", "nonlinearform", "restriction",
              "fespacehierarchy", "multigrid", "constraints",
              "transfer", "std_vectors",
              "tmop", "tmop_amr", "tmop_tools"]

    if add_cuda == '1':
        from setup_local import cudainc        
        include_dirs.append(cudainc)
    if add_libceed == '1':
        from setup_local import libceedinc        
        include_dirs.append(libceedinc)
    if add_suitesparse == '1':
        from setup_local import suitesparseinc        
        if suitesparseinc != "":
            include_dirs.append(suitesparseinc)
    if add_gslibs == '1':
        from setup_local import gslibsinc        
        include_dirs.append(gslibsinc)
        modules.append("gslib")

    sources = {name: [name + "_wrap.cxx"] for name in modules}
    proxy_names = {name: '_'+name for name in modules}

    tpl_include = []    
    for x in mfemstpl.split(' '):
        if x.startswith("-I"):
            tpl_include.append(x[2:])
    include_dirs.extend(tpl_include)
    
    extra_compile_args = [cxx11flag, '-DSWIG_TYPE_TABLE=PyMFEM']
    macros = [('TARGET_PY3', '1'), ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

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
                                  runtime_library_dirs = library_dirs,
                                  library_dirs = library_dirs,
                                  libraries = libraries,
                                  define_macros=macros)
                   for name in modules[1:]])

    return modules, ext_modules

def main():
    if 'clean' not in sys.argv:
        print('building serial version')                

    version = get_version()
    modules, ext_modules = get_extensions()

    setup (name = 'mfem_serial',
           version = version,
           author      = "S.Shiraiwa",
           description = """MFEM wrapper""",
           ext_modules = ext_modules,
           py_modules = modules, 
           )

if __name__ == '__main__':
    main()

