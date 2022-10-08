#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
print('building paralel version')

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
        import mpi4py
        mpi4pyinc = mpi4py.get_include()
    except ImportError:
        if 'clean' not in sys.argv:  raise        
        mpi4pyinc = ''
    
    try:
        from setup_local import (mfembuilddir, mfemincdir, mfemsrcdir, mfemlnkdir,
                                 mfemptpl,
                                 hypreinc, metisinc, hyprelib, metis5lib,
                                 cc_par, cxx_par,
                                 cxx11flag,                                 
                                 add_pumi, add_cuda, add_libceed, add_strumpack,
                                 add_suitesparse, add_gslibp)
        
        include_dirs = [mfembuilddir, mfemincdir, mfemsrcdir,
                        numpy.get_include(),
                        mpi4pyinc,
                        hypreinc, metisinc]
        library_dirs = [mfemlnkdir, hyprelib, metis5lib,]

    except ImportError:
        if 'clean' not in sys.argv:  raise
        cc_par = ''
        cxx_par = ''
        include_dirs = []
        library_dirs = []
        mfemptpl = ''        
        add_cuda = ''
        add_libceed = ''
        add_suitesparse = ''
        add_strumpack = ''
        add_pumi = ''
        add_gslibp = ''
        cxx11flag = ''

    libraries    = ['mfem', 'HYPRE', 'metis']
    
    ## remove current directory from path
    print("__file__", os.path.abspath(__file__))
    if '' in sys.path:
        sys.path.remove('')
    items = [x for x in sys.path if os.path.abspath(x) == os.path.dirname(os.path.abspath(__file__))]
    for x in items:
        sys.path.remove(x)
    print("sys path", sys.path)

    ## this forces to use compiler written in setup_local.py
    if cc_par != '': os.environ['CC'] = cc_par
    if cxx_par != '': os.environ['CXX'] = cxx_par

    modules= ["io_stream", "vtk", "sort_pairs", "datacollection",
              "globals", "mem_manager", "device", "hash", "stable3d",
              "cpointers", "symmat",
              "error", "array", "common_functions",
              "segment", "point", "hexahedron", "quadrilateral",
              "tetrahedron", "triangle", "wedge",
              "socketstream", "handle",
              "fe_base", "fe_fixed_order", "fe_h1", "fe_l2",
              "fe_nd", "fe_nurbs", "fe_pos", "fe_rt", "fe_ser", "doftrans",
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
              "fespacehierarchy", "multigrid", "constraints",
              "transfer", "dist_solver", "std_vectors", "auxiliary",
              "tmop", "tmop_amr", "tmop_tools"]

    if add_pumi == '1':
        from setup_local import puminc, pumilib
        modules.append("pumi")
        include_dirs.append(pumiinc)
        library_dirs.append(pumilib)

    if add_strumpack == '1':
        from setup_local import strumpackinc, strumpacklib
        modules.append("strumpack")
        if strumpackinc != "":
            include_dirs.append(strumpackinc)
        if strumpacklib != "":
            library_dirs.append(strumpacklib)

    if add_cuda == '1':
        from setup_local import cudainc
        include_dirs.append(cudainc)

    if add_libceed == '1':
        from setup_local import libceedinc
        include_dirs.append(libceedinc)

    if add_suitesparse == '1':
        from setup_local import suitesparseinc
        modules.append("schwarz")
        if suitesparseinc != "":
            include_dirs.append(suitesparseinc)

    if add_gslibp == '1':
        from setup_local import gslibpinc
        include_dirs.append(gslibpinc)
        modules.append("gslib")

    sources = {name: [name + "_wrap.cxx"] for name in modules}
    proxy_names = {name: '_'+name for name in modules}
        
    tpl_include = []    
    for x in mfemptpl.split(' '):
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
                            library_dirs = library_dirs,
                            runtime_library_dirs = library_dirs,
                            libraries = libraries,
                            define_macros=macros)
                   for name in modules[1:]])

    return modules, ext_modules

def main():
    if 'clean' not in sys.argv:
        print('building parallel version')                        

    version = get_version()
    modules, ext_modules = get_extensions()
    
    setup (name = 'mfem_parallel',
           version = version,
           author      = "S.Shiraiwa",
           description = """MFEM wrapper""",
           ext_modules = ext_modules,
           py_modules = modules, 
           )

if __name__ == '__main__':
    main()

