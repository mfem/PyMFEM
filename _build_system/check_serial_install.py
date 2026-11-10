"""
Check an installed release wheel and run example 1 in an isolated interpreter.

Run outside the checkout in a fresh venv without MPI dependencies:
python -I /path/to/check_serial_install.py /path/to/ex1.py /path/to/star.mesh
"""

import ctypes
import importlib.metadata
import importlib.util
import re
import runpy
import sys
from pathlib import Path


def loaded_libraries():
    """
    Report loaded libraries on the supported Linux and macOS wheel platforms.
    """
    if sys.platform == 'linux':
        return Path('/proc/self/maps').read_text()
    if sys.platform == 'darwin':
        dyld = ctypes.CDLL(None)
        dyld._dyld_image_count.restype = ctypes.c_uint32
        dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
        dyld._dyld_get_image_name.restype = ctypes.c_char_p
        return '\n'.join(dyld._dyld_get_image_name(i).decode()
                         for i in range(dyld._dyld_image_count()))
    raise RuntimeError('Serial wheel checks support Linux and macOS')


def main():
    """
    Check ownership, MPI absence, a real serial solve, and loaded libraries.
    """
    example, mesh = (str(Path(arg).resolve()) for arg in sys.argv[1:])
    assert sys.flags.isolated, 'Run this check with python -I'
    assert importlib.util.find_spec('mpi4py') is None, 'Use a venv without mpi4py'
    import mfem
    import mfem.ser

    assert Path(mfem.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert importlib.util.find_spec('mfem.par') is None
    assert importlib.util.find_spec('mfem._par') is None
    distribution = importlib.metadata.distribution('mfem')
    for entry in distribution.files:
        name = str(entry)
        assert not name.startswith(('mfem/_par/', 'mfem/external/par/')), name
        assert name != 'mfem/par.py', name
        assert distribution.locate_file(entry).is_file(), name
    print('Testing installed package:', mfem.__file__, flush=True)
    sys.argv = [example, '-m', mesh]
    runpy.run_path(example, run_name='__main__')
    assert 'mpi4py.MPI' not in sys.modules
    libraries = loaded_libraries().lower()
    assert not re.search(r'/lib(?:mpi|mpich|hypre)[^/\s]*', libraries), libraries
    print('Serial solve passed without loading MPI or HYPRE.')


if __name__ == '__main__':
    main()
