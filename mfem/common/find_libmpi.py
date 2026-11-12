"""
Report the MPI libraries loaded by mpi4py in this Python environment.

    Run: python -m mfem.common.find_libmpi [--check]

    Import mpi4py.MPI and return sorted, unique core MPI library paths.
    Report the runtime actually loaded, whether pip-installed or system MPI.
    (This does not inspect compiler-wrapper)
"""

import argparse
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

    return {'schema': 1, 'platform': platform.system(),
            'package': package.metadata['Name'].lower(), 'version': package.version,
            'provider': MPI.get_vendor()[0],
            'libraries': sorted({library_abi(path) for path in paths})}


def write_mpi_configuration(path, parallel=True):
    """
    Write parallel-build metadata, removing stale output for unchecked builds.
    """
    path = Path(path)
    path.unlink(missing_ok=True)
    config = mpi_configuration() if parallel else None
    if config is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(config, indent=2) + '\n', encoding='utf-8')
    return config


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
        if (not isinstance(expected, dict) or expected.get('schema') != 1 or
                expected.get('package') not in ('mpich', 'openmpi') or
                not all(isinstance(expected.get(key), str) and expected[key]
                        for key in ('platform', 'provider', 'version')) or
                not isinstance(expected.get('libraries'), list) or
                len(expected['libraries']) != 1 or
                not all(isinstance(name, str) and library_abi(name) == name
                        for name in expected['libraries'])):
            raise ValueError('Invalid MPI build configuration')
        if expected['platform'] != platform.system():
            raise RuntimeError(f'Expected platform {expected["platform"]}, '
                               f'detected {platform.system()}')
        actual = mpi_configuration()
        if actual is None:
            raise RuntimeError(f'Expected the {expected["package"]} Python package; none installed')
        for key in ('package', 'provider', 'libraries'):
            if actual[key] != expected[key]:
                raise RuntimeError(f'Expected MPI {key} {expected[key]!r}, detected {actual[key]!r}')
        return (f'MPI compatibility check passed: {actual["provider"]}, '
                f'{", ".join(actual["libraries"])}; '
                f'{actual["package"]} build version {expected["version"]}, '
                f'installed version {actual["version"]}.')
    except (ImportError, OSError, RuntimeError, ValueError) as error:
        raise ImportError(
            f'PyMFEM MPI compatibility check failed: {error}. '
            'Use an environment with the MPI Python package required by this wheel, '
            'or rebuild PyMFEM for your MPI. '
            'Diagnose with: python -m mfem.common.find_libmpi --check') from error


def main():
    """
    Print loaded paths or run the installed parallel wheel compatibility check.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Check the installed parallel build')
    args = parser.parse_args()
    try:
        print(check_mpi_runtime() if args.check else '\n'.join(loaded_mpi_libraries()))
    except (ImportError, OSError, RuntimeError) as error:
        parser.exit(1, f'MPI inspection failed: {error}\n')


if __name__ == '__main__':
    main()
