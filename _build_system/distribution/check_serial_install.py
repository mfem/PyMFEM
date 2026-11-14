"""
Check an installed release wheel and run example 1 in an isolated interpreter.

Run outside the checkout in a fresh venv without MPI dependencies:
python -I /path/to/check_serial_install.py /path/to/ex1.py /path/to/star.mesh
"""

import importlib.metadata
import importlib.util
import runpy
import sys
from pathlib import Path


def check_install_boundary():
    """
    Check the isolated base installation, file ownership and absent MPI payload.
    """
    assert sys.flags.isolated, 'Run this check with python -I'
    assert importlib.util.find_spec('mpi4py') is None, 'Use a venv without mpi4py'
    import mfem

    assert Path(mfem.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert importlib.util.find_spec('mfem.par') is None
    assert importlib.util.find_spec('mfem._par') is None
    distribution = importlib.metadata.distribution('mfem')
    for entry in distribution.files:
        name = str(entry)
        assert not name.startswith(('mfem/_par/', 'mfem/external/par/')), name
        assert name != 'mfem/par.py', name
        assert distribution.locate_file(entry).is_file(), name
    print('Serial installation boundary passed:', mfem.__file__, flush=True)


def run_serial_example(example, mesh):
    """
    Run example 1 with the installed serial bindings.
    """
    import mfem.ser  # noqa: F401 -- load the installed serial bindings

    sys.argv = [example, '-m', mesh]
    runpy.run_path(example, run_name='__main__')


def main():
    """
    Check the standalone installation, then run the serial example.
    """
    example, mesh = (str(Path(arg).resolve()) for arg in sys.argv[1:])
    check_install_boundary()
    run_serial_example(example, mesh)


if __name__ == '__main__':
    main()
