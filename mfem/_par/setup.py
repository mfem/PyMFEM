#!/usr/bin/env python

"""
setup.py file for SWIG example
"""


import sys
import os
from pathlib import Path

# this remove *.py in this directory to be imported from setuptools
# Github workflow (next import) fails without this, because it loads
# array.py in current directoy
sys.path.remove(os.path.abspath(os.path.dirname(sys.argv[0])))
from distutils.core import Extension, setup
from setuptools import setup

ddd = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
root = os.path.abspath(os.path.join(ddd, '..', '..'))
build_system_dir = os.path.join(root, '_build_system')
sys.path.insert(0, build_system_dir)
from build_generatedwrapperext import Build_StateO3
sys.path.pop(0)


def get_version():
    # read version number from __init__.py
    path = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(__file__))),
                        '..', '__init__.py')
    fid = open(path, 'r')
    lines = fid.readlines()
    fid.close()
    for x in lines:
        if x.strip().startswith('__version__'):
            version = eval(x.split('=')[-1].strip())
    return version


def get_extensions():
    # first load variables from PyMFEM_ROOT/setup_local.py
    sys.path.insert(0, root)
    try:
        from setup_local import (mfembuilddir, mfemincdir, mfemsrcdir, mfemlnkdir,
                                 mfemptpl, numpyinc, mpi4pyinc, mpiinc,
                                 hypreinc, metisinc, hyprelib, metis5lib,
                                 cc_par, cxx_par, cc_ser, cxx_ser,
                                 cxxstdflag, mfem_outside, build_miniapps,
                                 add_pumi, add_cuda, add_libceed, add_strumpack,
                                 add_suitesparse, add_gslibp, bdist_wheel_dir,
                                 enable_numba_swig_bridge)

        include_dirs = [mfembuilddir, mfemincdir, mfemsrcdir,
                        numpyinc,
                        mpi4pyinc,
                        hypreinc,
                        metisinc,]

        library_dirs = [mfemlnkdir,]
    except ImportError:
        if 'clean' not in sys.argv:
            raise
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
        cxxstdflag = '-std=c++17'
        mfem_outside = '0'
        build_miniapps = '0'
        mpiinc = ''
        enable_numba_swig_bridge = '0'

    libraries = ['mfem',]
    #if build_miniapps !=  '0':
    #    libraries.append("mfem-common")

    # remove current directory from path
    print("__file__", os.path.abspath(__file__))
    if '' in sys.path:
        sys.path.remove('')
    items = [x for x in sys.path if os.path.abspath(
        x) == os.path.dirname(os.path.abspath(__file__))]
    for x in items:
        sys.path.remove(x)

    # this forces to use compiler written in setup_local.py
    # if MPIINC (the location of MPI is given), serial compiler is used to build
    if mpiinc == '':
        if cc_par != '':
            os.environ['CC'] = cc_par
        if cxx_par != '':
            os.environ['CXX'] = cxx_par
    else:
        if cc_par != '':
            os.environ['CC'] = cc_ser
        if cxx_par != '':
            os.environ['CXX'] = cxx_ser

    modules = ["config",
               "io_stream", "vtk", "sort_pairs", "datacollection",
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
               "bilininteg", "intrules", "intrules_cut",
               "sparsemat", "densemat",
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
               "tmop", "tmop_amr", "tmop_tools", "qspace", "qfunction",
               "quadinterpolator", "quadinterpolator_face",
               "submesh", "transfermap", "staticcond", "sidredatacollection",
               "psubmesh", "ptransfermap", "enzyme",
               "attribute_sets", "arrays_by_name",
               "hyperbolic",  "complex_densemat",
               "bounds", "integrator", "ordering",
               "dpg", "particleset", "particlevector", "fe_pyramid",
               "multivector", "dgmassinv", "lor", "filteredsolver"]

    if mpiinc != '':
        include_dirs.append(mpiinc)

    if add_pumi == '1':
        from setup_local import puminc, pumilib
        modules.append("pumi")
        include_dirs.append(pumiinc)

    if add_strumpack == '1':
        from setup_local import strumpackinc, strumpacklib
        modules.append("strumpack")
        if strumpackinc != "":
            include_dirs.append(strumpackinc)

    if add_cuda == '1':
        from setup_local import cudainc
        include_dirs.append(cudainc)
    else:
        modules.append("navier_solver")

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

    if enable_numba_swig_bridge == '1' and 'clean' not in sys.argv:
        from nsb_rt.build_info import get_build_info

        bridge_info = get_build_info()
        artifacts = Path(ddd) / 'nsb_artifacts'
        director_root = artifacts / 'director' / 'cpp'
        helper = 'nsb_director_bindings'
        helper_sources = [
            str(Path(ddd) / (helper + '_wrap.cxx')),
            str(director_root / 'nsb_director.cpp'),
        ]
        modules.append(helper)
        sources[helper] = helper_sources
        proxy_names[helper] = '_' + helper

        state_name = 'nsb_state_bindings_ext'
        modules.append(state_name)
        sources[state_name] = [str(Path(ddd) / '_nsb_state_bindings_wrap.cxx'),
                               *bridge_info['sources']]
        proxy_names[state_name] = state_name
        missing = [path for name in (helper, state_name)
                   for path in sources[name] if not Path(path).is_file()]
        if missing:
            raise RuntimeError(
                'missing prepared PyMFEM NSB wrapper; run '
                'generate_nsb_pymfem_wrapper first: ' + missing[0])

    if enable_numba_swig_bridge == '1' and 'clean' not in sys.argv:
        mpi_name = '_nsb_mpi'
        modules.append(mpi_name)
        sources[mpi_name] = [str(artifacts / 'mpi' / 'cpp' / 'nsb_mpi.cpp')]
        proxy_names[mpi_name] = mpi_name
        include_dirs.extend([str(director_root), bridge_info['include_dir'],
                             bridge_info['cpp_dir']])

    tpl_include = []
    for x in mfemptpl.split(' '):
        if x.startswith("-I"):
            tpl_include.append(x[2:])
    include_dirs.extend(tpl_include)

    extra_compile_args = [cxxstdflag, '-DSWIG_TYPE_TABLE=PyMFEM']

    macros = [('TARGET_PY3', '1'),
              ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

    runtime_library_dirs = [
        x for x in library_dirs if x.find(bdist_wheel_dir) == -1]
    if mfem_outside == "0" and sys.platform  in ("linux", "linux2"):
        runtime_library_dirs.append("$ORIGIN/../external/par/lib")
        runtime_library_dirs.append("$ORIGIN/../external/lib")
    elif mfem_outside == "0" and sys.platform  == "darwin":
        runtime_library_dirs.append("@loader_path/../external/par/lib")
        runtime_library_dirs.append("@loader_path/../external/lib")
    else:
        runtime_library_dirs = library_dirs

    ext_modules = [Extension(proxy_names[modules[0]],
                             sources=sources[modules[0]],
                             extra_compile_args=extra_compile_args,
                             include_dirs=include_dirs,
                             library_dirs=library_dirs,
                             runtime_library_dirs=runtime_library_dirs,
                             libraries=libraries,
                             define_macros=macros), ]

    ext_modules.extend([Extension(proxy_names[name],
                                  sources=sources[name],
                                  extra_compile_args=extra_compile_args,
                                  include_dirs=include_dirs,
                                  library_dirs=library_dirs,
                                  runtime_library_dirs=runtime_library_dirs,
                                  libraries=libraries,
                                  define_macros=macros)
                        for name in modules[1:]])

    if enable_numba_swig_bridge == '1' and 'clean' not in sys.argv:
        state_extension = next(extension for extension in ext_modules
                               if extension.name == state_name)
        # Compile only the runtime implementation at -O3 on GCC/Clang.
        # The generated state wrapper retains the ordinary wrapper settings.
        state_extension.nsb_state_sources = tuple(bridge_info['sources'])
        state_extension.depends.extend([
            bridge_info['header'],
            str(Path(build_system_dir) / 'build_generatedwrapperext.py'),
        ])

    return modules, ext_modules


def main():
    if 'clean' not in sys.argv:
        print('building parallel version')

    version = get_version()
    modules, ext_modules = get_extensions()

    setup(name='mfem_parallel',
          version=version,
          author="S.Shiraiwa",
          description="""MFEM wrapper""",
          ext_modules=ext_modules,
          py_modules=modules,
          cmdclass={'build_ext': Build_StateO3},
          )


if __name__ == '__main__':
    main()
