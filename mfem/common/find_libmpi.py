"""
Inspect loaded MPI libraries and check compatibility with the wheel configuration.
"""

import ctypes
import json
import os
import platform
import re
from importlib import metadata
from pathlib import Path


def loaded_mpi_libraries():
    # load mpi4py and inspect what MPI library is dynamically linked to
    # Python executable.

    from mpi4py import MPI  # noqa: F401 -- import loads the MPI runtime

    def _linux_libraries():
        paths = []
        with Path('/proc/self/maps').open() as stream:
            for line in stream:
                fields = line.rstrip('\n').split(maxsplit=5)
                if len(fields) == 6 and fields[5].startswith('/'):
                    path = fields[5].removesuffix(' (deleted)')
                    paths.append(re.sub(r'\\([0-7]{3})',
                                    lambda match: chr(int(match[1], 8)), path))
        return paths

    def _macos_libraries():
        dyld = ctypes.CDLL(None)
        dyld._dyld_image_count.argtypes = []
        dyld._dyld_image_count.restype = ctypes.c_uint32
        dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
        dyld._dyld_get_image_name.restype = ctypes.c_char_p
        paths = []
        for index in range(dyld._dyld_image_count()):
            path = dyld._dyld_get_image_name(index)
            if path is not None:
                paths.append(os.fsdecode(path))
        return paths

    system = platform.system()
    if system not in ('Linux', 'Darwin'):
        raise RuntimeError(f'Unsupported MPI inspection platform: {system}')

    paths = _linux_libraries() if system == 'Linux' else _macos_libraries()
    libraries = sorted({path for path in paths if re.fullmatch(
        r'libmpi(?:_abi)?(?:\.so(?:\.\d+)*|(?:\.\d+)*\.dylib)', Path(path).name)})
    if not libraries:
        raise RuntimeError('No loaded libmpi library found after importing mpi4py.MPI')
    return libraries


def library_abi(path):
    """
    Return the core library name with its ABI major, excluding release suffixes.
    """
    name = Path(path).name
    match = re.fullmatch(r'(libmpi(?:_abi)?)\.so\.(\d+)(?:\.\d+)*', name)
    if match:
        return f'{match[1]}.so.{match[2]}'
    match = re.fullmatch(r'(libmpi(?:_abi)?)\.(\d+)(?:\.\d+)*\.dylib', name)
    if match:
        return f'{match[1]}.{match[2]}.dylib'
    raise RuntimeError(f'Cannot identify a versioned MPI ABI from {path}')


def mpi_configuration():
    """
    Describe loaded MPI and pip provider metadata.
    Return None when no pip MPI package is installed.
    """
    packages = []
    for name in ('mpich', 'openmpi'):
        try:
            packages.append(metadata.distribution(name))
        except metadata.PackageNotFoundError:
            pass
    if not packages:
        return None
    if len(packages) != 1:
        raise RuntimeError('Both mpich and openmpi Python packages are installed; keep one provider')
    package = packages[0]
    paths = loaded_mpi_libraries()
    if len(paths) != 1:
        raise RuntimeError(f'Expected one loaded core MPI library, found {paths}')
    from mpi4py import MPI

    return {'platform': platform.system(),
            'package': package.metadata['Name'].lower(), 'version': package.version,
            'provider': MPI.get_vendor()[0],
            'libraries': sorted({library_abi(path) for path in paths})}


def check_mpi_runtime(config_path=None):
    """
    Validate loaded pip MPI before native parallel imports; return a diagnostic.
    """
    path = (Path(config_path) if config_path is not None else
            Path(__file__).resolve().parents[1] / '_par' / 'mpi-runtime.json')
    try:
        try:
            text = path.read_text(encoding='utf-8')
        except FileNotFoundError:
            return 'No pip MPI build configuration installed; compatibility check skipped.'
        
        expected = json.loads(text)
        if (not isinstance(expected, dict) or expected.get('schema') != 2 or
                expected.get('package') not in ('mpich', 'openmpi') or
                not all(isinstance(expected.get(key), str) and expected[key]
                        for key in ('platform', 'provider', 'version')) or
                not isinstance(expected.get('libraries'), list) or
                not expected['libraries'] or
                not all(isinstance(name, str) and Path(name).name == name and
                        (name.endswith(('.so', '.dylib')) or '.so.' in name)
                        for name in expected['libraries'])):
            raise ValueError('Invalid MPI build configuration')
        
        if expected['platform'] != platform.system():
            raise RuntimeError(f'Expected platform {expected["platform"]}, '
                               f'detected {platform.system()}')
        actual = mpi_configuration()
        if actual is None:
            raise RuntimeError(f'Expected the {expected["package"]} Python package; none installed')
        for key in ('package', 'provider'):
            if actual[key] != expected[key]:
                raise RuntimeError(f'Expected MPI {key} {expected[key]!r}, detected {actual[key]!r}')
        if not set(actual['libraries']).issubset(expected['libraries']):
            raise RuntimeError(f'Expected MPI libraries supplied by {expected["package"]}: '
                               f'{expected["libraries"]!r}, detected {actual["libraries"]!r}')
        return (f'MPI compatibility check passed: {actual["provider"]}, '
                f'{", ".join(actual["libraries"])}; '
                f'{actual["package"]} build version {expected["version"]}, '
                f'installed version {actual["version"]}.')
    
    except (ImportError, OSError, RuntimeError, ValueError) as error:
        raise ImportError(
            f'PyMFEM MPI compatibility check failed: {error}. '
            'Use an environment with the MPI Python package required by this wheel, '
            'or rebuild PyMFEM for your MPI. '
            'Diagnose by calling mfem.common.find_libmpi.check_mpi_runtime().') from error
