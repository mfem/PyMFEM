"""
Exercise a composed base/add-on installation outside the source checkout.

Run each mode with python -I in a separate process. The parallel mode requires
mpiexec -n 2 and runs example 1p followed by an explicit GSLIB operation.
"""

import argparse
import importlib
import runpy
import sys
from pathlib import Path

# Isolated Python omits the script directory from sys.path. Load the existing
# checker explicitly without adding the checkout to the import search path.
loaded_libraries = runpy.run_path(
    str(Path(__file__).with_name('check_serial_install.py')))['loaded_libraries']


def check_runtime():
    """
    Require the installed MPICH runtime and reject build-time MPI paths.
    """
    from mpi4py import MPI

    assert 'MPICH' in MPI.Get_library_version(), MPI.Get_library_version()
    assert MPI.COMM_WORLD.size == 2, 'Run the parallel check with two ranks'
    if sys.platform == 'linux':
        paths = {line.split()[-1] for line in loaded_libraries().splitlines()
                 if '/libmpi' in line}
        assert paths, 'No MPI runtime is loaded'
        assert all(Path(path).resolve().is_relative_to(Path(sys.prefix).resolve())
                   for path in paths), paths
        assert all('libmpi.so.12' in Path(path).name for path in paths), paths
    return MPI


def run_example(example, mesh):
    """
    Run an existing solve using an explicit absolute mesh path.
    """
    sys.argv = [str(example), '-m', str(mesh)]
    runpy.run_path(str(example), run_name='__main__')


def check_gslib(MPI):
    """
    Exercise point finding explicitly so an optional import cannot hide failure.
    """
    import mfem.par as mfem
    gslib = importlib.import_module('mfem._par.gslib')

    mesh = mfem.Mesh.MakeCartesian2D(2, 2, mfem.Element.QUADRILATERAL, True, 1., 1.)
    mesh.SetCurvature(1)
    parallel_mesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    finder = gslib.FindPointsGSLIB(MPI.COMM_WORLD)
    finder.Setup(parallel_mesh)
    points = mfem.Vector([0.25, 0.25])
    finder.FindPoints(points)
    codes = finder.GetCode()
    assert codes.Size() == 1 and codes[0] < 2, codes.ToList()
    finder.FreeData()
    print(f'Rank {MPI.COMM_WORLD.rank}: parallel GSLIB passed', flush=True)


def main():
    """
    Check serial composition, both entry paths, and both import guard orders.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['serial', 'parallel', 'direct',
                                        'serial-first', 'parallel-first'])
    parser.add_argument('--example', type=Path)
    parser.add_argument('--mesh', type=Path)
    args = parser.parse_args()
    if args.mode in ('serial', 'parallel') and (args.example is None or args.mesh is None):
        parser.error('--example and --mesh are required for serial/parallel solves')
    assert sys.flags.isolated, 'Run this check with python -I'
    import mfem
    assert Path(mfem.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    if args.mode in ('serial-first', 'parallel-first'):
        first, second = ('ser', 'par') if args.mode == 'serial-first' else ('par', 'ser')
        importlib.import_module('mfem.' + first)
        try:
            importlib.import_module('mfem.' + second)
        except ImportError as error:
            assert 'mode is already loaded' in str(error), error
        else:
            raise AssertionError('Serial/parallel import guard did not fire')
    elif args.mode == 'direct':
        import mfem._par
        assert 'mpi4py.MPI' in sys.modules
        # Exercise a direct extension import after the package initialized MPI.
        from mfem._par.mesh import Mesh
        assert Mesh is not None
    elif args.mode == 'parallel':
        # mfem.par must initialize MPI before its own extension imports.
        import mfem.par
        MPI = check_runtime()
        run_example(args.example.resolve(), args.mesh.resolve())
        check_gslib(MPI)
        check_runtime()
    else:
        import mfem.ser
        run_example(args.example.resolve(), args.mesh.resolve())
        assert 'mpi4py.MPI' not in sys.modules
        assert '/libmpi' not in loaded_libraries().lower()
    print(f'Composed installation check passed: {args.mode}', flush=True)


if __name__ == '__main__':
    main()
