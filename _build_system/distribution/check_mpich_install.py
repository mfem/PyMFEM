"""
Run a two-rank parallel solve and GSLIB operation with the installed wheels.

Stage ex1p.py and star.mesh beside this script, then run:
mpiexec -n 2 python -I /path/to/check_mpich_install.py
"""

import importlib
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def check_runtime():
    """
    Require the installed MPICH runtime and reject build-time MPI paths.
    """
    from mpi4py import MPI

    from mfem.common.find_libmpi import loaded_mpi_libraries

    assert 'MPICH' in MPI.Get_library_version(), MPI.Get_library_version()
    assert MPI.COMM_WORLD.size == 2, 'Run the parallel check with two ranks'
    paths = loaded_mpi_libraries()
    assert all(Path(path).resolve().is_relative_to(Path(sys.prefix).resolve())
               for path in paths), paths
    return MPI


def run_example():
    """
    Run an existing solve using an explicit absolute mesh path.
    """
    sys.argv = [str(HERE / 'ex1p.py'), '-m', str(HERE / 'star.mesh')]
    runpy.run_path(sys.argv[0], run_name='__main__')


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
    Check the parallel entry point, runtime, solve and explicit GSLIB operation.
    """
    assert sys.flags.isolated, 'Run this check with python -I'
    import mfem.par  # noqa: F401
    MPI = check_runtime()
    run_example()
    check_gslib(MPI)
    check_runtime()


if __name__ == '__main__':
    main()
