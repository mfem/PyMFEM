# ----------------------------------------------------------------------------------------
# Routines for PyMFEM Wrapper Generation/Compile
# ----------------------------------------------------------------------------------------
import sys
import os
import re
import subprocess

__all__ = ["write_setup_local", "generate_wrapper", "generate_nsb_artifacts",
           "generate_nsb_pymfem_wrapper",
           "clean_wrapper", "make_mfem_wrapper"]

from build_utils import *
from build_consts import *
import build_globals as bglb


def write_setup_local():
    '''
    create setup_local.py. parameters written here will be read
    by setup.py in mfem._ser and mfem._par
    '''
    mfemser = bglb.mfems_prefix
    mfempar = bglb.mfemp_prefix

    hyprelibpath = os.path.dirname(
        find_libpath_from_prefix('HYPRE', bglb.hypre_prefix))
    metislibpath = os.path.dirname(
        find_libpath_from_prefix('metis', bglb.metis_prefix))

    mfems_tpl = read_mfem_tplflags(bglb.mfems_prefix)
    mfemp_tpl = read_mfem_tplflags(
        bglb.mfemp_prefix) if bglb.build_parallel else ''

    print(mfems_tpl, mfemp_tpl)

    params = {'cxx_ser': bglb.cxx_command,
              'cc_ser': bglb.cc_command,
              'cxx_par': bglb.mpicxx_command,
              'cc_par': bglb.mpicc_command,
              'whole_archive': '--whole-archive',
              'no_whole_archive': '--no-whole-archive',
              'nocompactunwind': '',
              'swigflag': '-Wall -c++ -python -fastproxy -olddefs -keyword',
              'hypreinc': os.path.join(bglb.hypre_prefix, 'include'),
              'hyprelib': hyprelibpath,
              'metisinc': os.path.join(bglb.metis_prefix, 'include'),
              'metis5lib': metislibpath,
              'numpyinc': get_numpy_inc(),
              'mpi4pyinc': '',
              'mpiinc':bglb.mpiinc,
              'mfem_outside': '1' if  bglb.mfem_outside else '0',
              'mfembuilddir': os.path.join(mfempar, 'include'),
              'mfemincdir': os.path.join(mfempar, 'include', 'mfem'),
              'mfemlnkdir': os.path.join(mfempar, 'lib'),
              'mfemserbuilddir': os.path.join(mfemser, 'include'),
              'mfemserincdir': os.path.join(mfemser, 'include', 'mfem'),
              'mfemserlnkdir': os.path.join(mfemser, 'lib'),
              'mfemsrcdir': os.path.join(bglb.mfem_source),
              'mfemstpl': mfems_tpl,
              'mfemptpl': mfemp_tpl,
              'add_pumi': '',
              'add_strumpack': '',
              'add_cuda': '',
              'add_libceed': '',
              'add_suitesparse': '',
              'add_gslib': '',
              'add_gslibp': '',
              'add_gslibs': '',
              'libceedinc': os.path.join(bglb.libceed_prefix, 'include'),
              'gslibsinc': os.path.join(bglb.gslibs_prefix, 'include'),
              'gslibpinc': os.path.join(bglb.gslibp_prefix, 'include'),
              'cxxstdflag': bglb.cxxstd_flag,
              'build_mfem': '1' if bglb.build_mfem else '0',
              'build_miniapps': '1' if bglb.mfem_miniapps else '0',
              'enable_numba_swig_bridge': ('1' if bglb.enable_numba_swig_bridge
                                           else '0'),
              'bdist_wheel_dir': bglb.bdist_wheel_dir,
              }

    if bglb.build_parallel:
        params['mpi4pyinc'] = get_mpi4py_inc()

    def add_extra(xxx, inc_sub=None):
        params['add_' + xxx] = '1'
        ex_prefix = getattr(bglb, xxx + '_prefix')
        if inc_sub is None:
            params[xxx +
                   'inc'] = os.path.join(ex_prefix, 'include')
        else:
            params[xxx +
                   'inc'] = os.path.join(ex_prefix, 'include', inc_sub)

        params[xxx + 'lib'] = os.path.join(ex_prefix, 'lib')

    if bglb.enable_pumi:
        add_extra('pumi')
    if bglb.enable_strumpack:
        add_extra('strumpack')
    if bglb.enable_cuda:
        add_extra('cuda')
    if bglb.enable_libceed:
        add_extra('libceed')
    if bglb.enable_suitesparse:
        add_extra('suitesparse', inc_sub='suitesparse')
    if bglb.enable_gslib:
        add_extra('gslibs')
    if bglb.enable_gslib:
        add_extra('gslibp')

    pwd = chdir(rootdir)

    fid = open('setup_local.py', 'w')
    fid.write("#  setup_local.py \n")
    fid.write("#  generated from setup.py\n")
    fid.write("#  do not edit this directly\n")

    for key, value in params.items():
        text = key.lower() + ' = "' + value + '"'
        fid.write(text + "\n")
    fid.close()

    os.chdir(pwd)


def generate_wrapper(do_parallel):
    '''
    run swig.
    '''
    # this should work as far as we are in the same directory ?
    from multiprocessing import Pool
    import build_globals as bglb

    if bglb.dry_run or bglb.verbose:
        print("generating SWIG wrapper")
        print("using MFEM source", os.path.abspath(bglb.mfem_source))
    if not os.path.exists(os.path.abspath(bglb.mfem_source)):
        assert False, "MFEM source directory. Use --mfem-source=<path>"

    def ifiles():
        ifiles = os.listdir()
        ifiles = [x for x in ifiles if x.endswith('.i')]
        ifiles = [x for x in ifiles if not x.startswith('#')]
        ifiles = [x for x in ifiles if not x.startswith('.')]
        return ifiles

    def check_new(ifile):
        wfile = ifile[:-2] + '_wrap.cxx'
        if not os.path.exists(wfile):
            return True
        return os.path.getmtime(ifile) > os.path.getmtime(wfile)

    def update_integrator_exts():
        pwd = chdir(os.path.join(rootdir, 'mfem', 'common'))
        command1 = [sys.executable, "generate_lininteg_ext.py"]
        command2 = [sys.executable, "generate_bilininteg_ext.py"]
        make_call(command1)
        make_call(command2)
        os.chdir(pwd)

    def update_header_exists(mfem_source):
        print("updating the list of existing headers")
        list_of_headers = []
        L = len(mfem_source.split(os.sep))
        for (dirpath, dirnames, filenames) in os.walk(mfem_source):
            for filename in filenames:
                if filename.endswith('.hpp'):
                    dirs = dirpath.split(os.sep)[L:]
                    dirs.append(filename[:-4])
                    tmp = '_'.join(dirs)
                    xx = re.split('_|-', tmp)
                    new_name = 'FILE_EXISTS_'+'_'.join([x.upper() for x in xx])
                    if new_name not in list_of_headers:
                        list_of_headers.append(new_name)

        pwd = chdir(os.path.join(rootdir, 'mfem', 'common'))
        fid = open('existing_mfem_headers.i', 'w')
        for x in list_of_headers:
            fid.write("#define " + x + "\n")
        fid.close()
        os.chdir(pwd)

    mfemser = bglb.mfems_prefix
    mfempar = bglb.mfemp_prefix

    update_header_exists(bglb.mfem_source)

    swigflag = '-Wall -c++ -python -std=c++17 -fastproxy -olddefs -keyword'.split(' ')
    bridgeflag = []
    bridge_include_dir = None
    if bglb.enable_numba_swig_bridge:
        # The backend installs this build requirement only when the matching
        # PEP 517 setting is enabled. Its include directory provides the
        # annotation macros used by mfem/common/bridges/*.i.
        from nsb_rt.build_info import get_build_info
        bridge_build_info = get_build_info()
        bridge_include_dir = bridge_build_info['include_dir']
        bridgeflag = ['-DNUMBA_SWIG_BRIDGE', '-I' + bridge_include_dir]

        def committed_serial_bridge_flags():
            """Use the reviewed serial NSB artifacts committed in the source tree."""
            from pathlib import Path
            import json

            artifacts = Path(rootdir) / "mfem" / "_ser" / "nsb_artifacts"
            manifest_path = artifacts / "manifest.json"
            if not manifest_path.is_file():
                raise RuntimeError("missing committed serial NSB artifact manifest: " + str(manifest_path))
            manifest = json.loads(manifest_path.read_text())
            if manifest.get("schema_version") != 2 or manifest.get("artifact_schema_version") != 1:
                raise RuntimeError("unsupported committed serial NSB artifact schema")
            if manifest.get("state_backend") != "mfem":
                raise RuntimeError("serial PyMFEM NSB artifacts must use the 'mfem' backend")
            return ["-I" + str(artifacts), "-DNUMBA_SWIG_BRIDGE_GENERATED",
                    "-DMFEM_NUMBA_SWIG_BRIDGE_SER"]
        serial_bridge_flags = committed_serial_bridge_flags()
        serial_bridge_modules = (
            "vector", "densemat", "array", "doftrans", "fespace", "fe_base",
            "eltrans", "intrules", "coefficient", "lininteg", "nonlinearform",
            "hyperbolic", "operators",
        )
        serial_bridge_interfaces = {module_name + ".i"
                                    for module_name in serial_bridge_modules}
    else:
        serial_bridge_flags = []
        serial_bridge_interfaces = set()

    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))

    serflag = ['-I' + os.path.join(mfemser, 'include'),
               '-I' + os.path.join(mfemser, 'include', 'mfem'),
               '-I' + os.path.abspath(bglb.mfem_source)]
    if bglb.enable_suitesparse:
        serflag.append('-I' + os.path.join(bglb.suitesparse_prefix,
                                           'include', 'suitesparse'))

    for filename in ['lininteg.i', 'bilininteg.i']:
        interface_bridge = (serial_bridge_flags
                            if filename in serial_bridge_interfaces else [])
        command = ([swig_command] + swigflag + bridgeflag + interface_bridge
                   + serflag + [filename])
        make_call(command)
    update_integrator_exts()

    commands = []
    for filename in ifiles():
        if not check_new(filename):
            continue
        interface_bridge = (serial_bridge_flags
                            if filename in serial_bridge_interfaces else [])
        command = [swig_command] + swigflag + bridgeflag + interface_bridge + serflag + [filename]
        commands.append(command)

    mp_pool = Pool(max((cpu_count() - 1, 1)))
    with mp_pool:
        mp_pool.map(subprocess.run, commands)

    if not do_parallel:
        os.chdir(pwd)
        return

    chdir(os.path.join(rootdir, 'mfem', '_par'))

    parflag = ['-I' + os.path.join(mfempar, 'include'),
               '-I' + os.path.join(mfempar, 'include', 'mfem'),
               '-I' + os.path.abspath(bglb.mfem_source),
               '-I' + os.path.join(bglb.hypre_prefix, 'include'),
               '-I' + os.path.join(bglb.metis_prefix, 'include'),
               '-I' + get_mpi4py_inc()]

    if bglb.enable_pumi:
        parflag.append('-I' + os.path.join(bglb.pumi_prefix, 'include'))
    if bglb.enable_strumpack:
        parflag.append('-I' + os.path.join(bglb.strumpack_prefix, 'include'))
    if bglb.enable_suitesparse:
        parflag.append('-I' + os.path.join(bglb.suitesparse_prefix,
                                           'include', 'suitesparse'))
    commands = []
    for filename in ifiles():
        if filename == 'strumpack.i' and not bglb.enable_strumpack:
            continue
        if not check_new(filename):
            continue
        command = [swig_command] + swigflag + bridgeflag + parflag + [filename]
        commands.append(command)

    mp_pool = Pool(max((cpu_count() - 1, 1)))
    with mp_pool:
        mp_pool.map(subprocess.run, commands)

    os.chdir(pwd)


def generate_nsb_pymfem_wrapper():
    """Prepare serial PyMFEM NSB wrappers from committed artifacts."""
    if not bglb.enable_numba_swig_bridge:
        return

    import json
    import shutil
    from pathlib import Path
    from nsb_rt.build_info import get_build_info

    root = Path(rootdir)
    package_dir = root / 'mfem' / '_ser'
    artifacts = package_dir / 'nsb_artifacts'
    manifest_path = artifacts / 'manifest.json'
    if not manifest_path.is_file():
        raise RuntimeError('missing committed serial NSB artifact manifest: ' + str(manifest_path))
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema_version') != 2 or manifest.get('artifact_schema_version') != 1:
        raise RuntimeError('unsupported committed serial NSB artifact schema')
    if manifest.get('state_backend') != 'mfem':
        raise RuntimeError("serial PyMFEM NSB artifacts must use the 'mfem' backend")

    bridge_info = get_build_info()
    required_runtime = manifest.get('required_runtime_version')
    runtime_version = bridge_info.get('version')
    try:
        required_parts = tuple(int(part) for part in required_runtime.split('.'))
        runtime_parts = tuple(int(part) for part in runtime_version.split('.'))
    except (AttributeError, TypeError, ValueError) as error:
        raise RuntimeError('invalid PyMFEM NSB runtime version metadata') from error
    if runtime_parts < required_parts:
        raise RuntimeError(
            'serial PyMFEM NSB artifacts require runtime >= '
            f'{required_runtime}, but build uses {runtime_version}'
        )
    if manifest.get('state_abi') != bridge_info['abi_version']:
        raise RuntimeError('serial PyMFEM NSB artifacts require a different state ABI')

    # Check the complete committed set before generating anything.  This keeps
    # a partial or stale artifact commit from producing a seemingly valid build.
    generated_files = manifest.get('generated_files')
    if not isinstance(generated_files, list) or 'manifest.json' not in generated_files:
        raise RuntimeError('committed serial NSB manifest has no generated file inventory')
    missing_generated = [artifacts / relative for relative in generated_files
                         if not (artifacts / relative).is_file()]
    if missing_generated:
        raise RuntimeError('missing committed serial NSB artifact: ' + str(missing_generated[0]))

    inventory_path = artifacts / manifest.get('inventory', '')
    if not inventory_path.is_file():
        raise RuntimeError('missing committed serial NSB inventory: ' + str(inventory_path))
    inventory = json.loads(inventory_path.read_text())
    sources = inventory.get('sources', [])
    expected_sources = {
        (entry.get('module'), entry.get('swig_module'), entry.get('interface'))
        for entry in sources
    }
    if len(expected_sources) != len(sources) or not sources:
        raise RuntimeError('invalid committed serial NSB source inventory')
    manifest_extensions = manifest.get('extensions', {})
    inventory_modules = {entry[0] for entry in expected_sources}
    if set(manifest_extensions) != inventory_modules:
        raise RuntimeError('serial NSB manifest and inventory source sets differ')
    for module, swig_module, interface in expected_sources:
        interface_path = root / interface
        extension = manifest_extensions[module]
        helper_path = artifacts / extension.get('helper_fragment', '')
        if not interface_path.is_file() or not helper_path.is_file():
            raise RuntimeError('serial NSB source or helper fragment is missing for ' + module)
        if extension.get('interface') != interface:
            raise RuntimeError('serial NSB manifest interface mismatch for ' + module)
        if extension.get('helper_fragment') not in generated_files:
            raise RuntimeError('serial NSB helper fragment is not in generated file inventory for ' + module)

    for descriptor_key in ('ordinary_descriptor', 'director_descriptor'):
        descriptor_path = artifacts / manifest.get(descriptor_key, '')
        if not descriptor_path.is_file():
            raise RuntimeError('missing committed serial NSB descriptor: ' + str(descriptor_path))
        descriptor = json.loads(descriptor_path.read_text())
        if descriptor.get('bridge_id') != manifest.get('bridge_id'):
            raise RuntimeError('serial NSB descriptor bridge ID mismatch: ' + str(descriptor_path))
        descriptor_sources = descriptor.get('source_inventory', {}).get('sources')
        if descriptor_sources != sources:
            raise RuntimeError('serial NSB descriptor source inventory mismatch: ' + str(descriptor_path))

    for relative in manifest.get('python_files', []):
        source = artifacts / relative
        destination = root / Path(relative).relative_to('python')
        if not source.is_file():
            raise RuntimeError('missing committed serial NSB artifact: ' + str(source))
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    director = manifest['director']
    director_root = artifacts / Path(director['sources'][0]).parent
    director_wrapper = package_dir / 'nsb_director_bindings_wrap.cxx'
    director_includes = [
        str(director_root), str(package_dir), bridge_info['include_dir'],
        bridge_info['cpp_dir'], os.path.join(bglb.mfems_prefix, 'include'),
        os.path.join(bglb.mfems_prefix, 'include', 'mfem'),
        os.path.abspath(bglb.mfem_source),
    ]
    make_call([
        swig_command, '-c++', '-python',
        *['-I' + path for path in director_includes],
        '-outdir', str(package_dir), '-o', str(director_wrapper),
        str(artifacts / director['swig_interface']),
    ])

    state_wrapper = package_dir / '_nsb_state_bindings_wrap.cxx'
    state_proxy_dir = root / 'mfem'
    state_includes = [bridge_info['include_dir'], bridge_info['cpp_dir']]
    make_call([
        swig_command, '-c++', '-python', '-module', '_nsb_state_bindings',
        '-interface', '_nsb_state_bindings_ext',
        *['-I' + path for path in state_includes],
        '-outdir', str(state_proxy_dir), '-o', str(state_wrapper),
        bridge_info['state_interface'],
    ])


def generate_nsb_artifacts(check=False):
    """Regenerate or check PyMFEM's committed NSB artifact directory.

    This is an explicit developer operation.  Ordinary builds call
    ``generate_nsb_pymfem_wrapper`` and consume the committed artifacts; they
    do not import the private generator package.
    """
    from pathlib import Path

    config = Path(rootdir) / 'mfem' / '_ser' / 'bridge.toml'
    if not config.is_file():
        raise RuntimeError('missing PyMFEM NSB generator configuration: ' + str(config))
    command = [sys.executable, '-m', 'numba_swig_bridge.generator',
               '--config', str(config)]
    if check:
        command.append('--check')
    make_call(command, force_verbose=True)


def clean_wrapper():
    from pathlib import Path

    # serial
    pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.cxx')]

    print(os.getcwd(), wfiles)
    remove_files(wfiles)

    wfiles = [x for x in os.listdir() if x.endswith('_wrap.h')]
    remove_files(wfiles)

    wfiles = [x for x in os.listdir() if x.endswith('.py')]
    wfiles.remove("__init__.py")
    wfiles.remove("setup.py")
    wfiles.remove("tmop_modules.py")
    if "bridge_registration.py" in wfiles:
        wfiles.remove("bridge_registration.py")
    remove_files(wfiles)

    ifiles = [x for x in os.listdir() if x.endswith('.i')]
    for x in ifiles:
        Path(x).touch()

    # parallel
    chdir(os.path.join(rootdir, 'mfem', '_par'))
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.cxx')]

    remove_files(wfiles)
    wfiles = [x for x in os.listdir() if x.endswith('_wrap.h')]
    remove_files(wfiles)

    wfiles = [x for x in os.listdir() if x.endswith('.py')]
    wfiles.remove("__init__.py")
    wfiles.remove("setup.py")
    wfiles.remove("tmop_modules.py")
    if "bridge_registration.py" in wfiles:
        wfiles.remove("bridge_registration.py")
    remove_files(wfiles)

    ifiles = [x for x in os.listdir() if x.endswith('.i')]
    for x in ifiles:
        Path(x).touch()

    chdir(pwd)


def make_mfem_wrapper(serial=True):
    '''
    compile PyMFEM wrapper code
    '''
    import build_globals as bglb

    if bglb.dry_run or bglb.verbose:
        print("compiling wrapper code, serial=" + str(serial))
    if not os.path.exists(os.path.abspath(bglb.mfem_source)):
        assert False, "MFEM source directory. Use --mfem-source=<path>"

    record_mfem_sha(bglb.mfem_source)

    write_setup_local()

    if serial:
        generate_nsb_pymfem_wrapper()
        pwd = chdir(os.path.join(rootdir, 'mfem', '_ser'))
    else:
        pwd = chdir(os.path.join(rootdir, 'mfem', '_par'))

    python = sys.executable


    command = [python, 'setup.py', 'build_ext', '--inplace',
               '--parallel',  str(cpu_count())]

    state_package_link = None
    if serial and bglb.enable_numba_swig_bridge:
        # The serial setup script runs in mfem/_ser, but the qualified state
        # extension belongs in the parent mfem package.  A temporary package
        # link lets setuptools resolve that qualified extension correctly.
        state_package_link = os.path.join(os.getcwd(), 'mfem')
        if not os.path.lexists(state_package_link):
            os.symlink('..', state_package_link)
    try:
        make_call(command, force_verbose=True)
    finally:
        if state_package_link and os.path.islink(state_package_link):
            os.unlink(state_package_link)

    os.chdir(pwd)
