#!/usr/bin/env python

"""
Serial version setup file
"""


import sys
import os
from pathlib import Path

# this remove *.py in this directory to be imported from setuptools
# Github workflow (next import) fails without this, because it loads
# array.py in current directoy
sys.path.remove(os.path.abspath(os.path.dirname(sys.argv[0])))
from distutils.core import Extension, setup

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
        from setup_local import (mfemserbuilddir, mfemserincdir, mfemsrcdir,
                                 mfemserlnkdir, mfemstpl, numpyinc,
                                 cc_ser, cxx_ser,
                                 cxxstdflag, mfem_outside, build_miniapps,
                                 add_cuda, add_libceed, add_suitesparse, add_gslibs,
                                 enable_numba_swig_bridge,
                                 bdist_wheel_dir)

        include_dirs = [mfemserbuilddir, mfemserincdir, mfemsrcdir, numpyinc,]
        library_dirs = [mfemserlnkdir, ]

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
        cxxstdflag = '-std=c++17'
        mfem_outside = '0'
        build_miniapps = '0'
        enable_numba_swig_bridge = '0'


    libraries = ['mfem']
    #if build_miniapps != '0':
    #    libraries.append("mfem-common")

    # remove current directory from path
    # print("__file__", os.path.abspath(__file__))
    if '' in sys.path:
        sys.path.remove('')
    items = [x for x in sys.path if os.path.abspath(
        x) == os.path.dirname(os.path.abspath(__file__))]
    for x in items:
        sys.path.remove(x)
    # print("sys path", sys.path)

    # this forces to use compiler written in setup_local.py
    if cc_ser != '':
        os.environ['CC'] = cc_ser
    if cxx_ser != '':
        os.environ['CXX'] = cxx_ser

    modules = ["config",
               "io_stream", "vtk", "sort_pairs", "datacollection",
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
               "bilininteg", "intrules", "intrules_cut",
               "sparsemat", "densemat",
               "solvers", "estimators", "mesh_operators", "ode",
               "sparsesmoothers",
               "matrix", "operators", "ncmesh", "eltrans", "geom",
               "nonlininteg", "nonlinearform", "restriction",
               "fespacehierarchy", "multigrid", "constraints",
               "transfer", "std_vectors",
               "tmop", "tmop_amr", "tmop_tools", "qspace", "qfunction",
               "quadinterpolator", "quadinterpolator_face",
               "submesh", "transfermap", "staticcond",
               "sidredatacollection", "enzyme",
               "attribute_sets", "arrays_by_name",
               "hyperbolic", "complex_densemat", 
               "bounds", "integrator", "ordering", 
               "dpg", "particleset", "particlevector", "fe_pyramid",
               "multivector", "dgmassinv", "lor", "filteredsolver"]

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

        state_name = 'mfem._nsb_state_bindings_ext'
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

    tpl_include = []
    for x in mfemstpl.split(' '):
        if x.startswith("-I"):
            if x.find("MacOS") != -1 and x.find(".sdk") != -1:
                continue
            tpl_include.append(x[2:])
    include_dirs.extend(tpl_include)

    if enable_numba_swig_bridge == '1' and 'clean' not in sys.argv:
        include_dirs.extend([str(director_root), bridge_info['include_dir'],
                             bridge_info['cpp_dir']])

    extra_compile_args = [cxxstdflag, '-DSWIG_TYPE_TABLE=PyMFEM']
    macros = [('TARGET_PY3', '1'),
              ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

    runtime_library_dirs = [
        x for x in library_dirs if x.find(bdist_wheel_dir) == -1]
    if mfem_outside == "0" and sys.platform in ("linux", "linux2"):
        runtime_library_dirs.append("$ORIGIN/../external/ser/lib")
    elif mfem_outside == "0" and sys.platform == "darwin":
        runtime_library_dirs.append("@loader_path/../external/ser/lib")
    else:
        runtime_library_dirs = library_dirs

    ext_modules = [Extension(proxy_names[modules[0]],
                             sources=sources[modules[0]],
                             extra_compile_args=extra_compile_args,
                             extra_link_args=[],
                             include_dirs=include_dirs,
                             library_dirs=library_dirs,
                             runtime_library_dirs=runtime_library_dirs,
                             libraries=libraries,
                             define_macros=macros), ]

    ext_modules.extend([Extension(proxy_names[name],
                                  sources=sources[name],
                                  extra_compile_args=extra_compile_args,
                                  extra_link_args=[],
                                  include_dirs=include_dirs,
                                  runtime_library_dirs=runtime_library_dirs,
                                  library_dirs=library_dirs,
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
        print('building serial version')

    version = get_version()
    modules, ext_modules = get_extensions()
    python_modules = [name for name in modules if '.' not in name]

    setup(name='mfem_serial',
          version=version,
          author="S.Shiraiwa",
          description="""MFEM wrapper""",
          ext_modules=ext_modules,
          py_modules=python_modules,
          cmdclass={'build_ext': Build_StateO3},
          )


if __name__ == '__main__':
    main()
