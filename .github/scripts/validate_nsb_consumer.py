"""Run NSB fixtures against a wheel in a separate consumer without build tools."""

import argparse
import json
import os
from pathlib import Path
import shutil
import shlex
import signal
import subprocess
import sys
import venv


ROOT = Path(__file__).resolve().parents[2]


def run(command, *, cwd, env, timeout=600):
    command = [str(value) for value in command]
    print("$ " + " ".join(command), flush=True)
    with subprocess.Popen(command, cwd=cwd, env=env, start_new_session=True) as child:
        try:
            status = child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
            raise
        if status:
            # Kill any surviving MPI workers before reporting the failure.
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            raise subprocess.CalledProcessError(status, command)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--backend', choices=('serial', 'openmpi', 'mpich'), required=True
    )
    parser.add_argument('--runtime-version', required=True)
    parser.add_argument('--wheel-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    wheels = args.wheel_dir.resolve()
    if work == ROOT or ROOT in work.parents:
        parser.error('consumer directory must be outside the checkout')
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        parser.error('consumer directory must be empty')
    pymfem_wheels = list(wheels.glob('mfem-*.whl'))
    if len(pymfem_wheels) != 1:
        parser.error('expected exactly one PyMFEM wheel')
    mpi = args.backend != 'serial'
    env = dict(os.environ)
    for key in ('PYTHONPATH', 'PYTHONHOME', 'PIP_TARGET', 'PIP_PREFIX', 'PIP_USER'):
        env.pop(key, None)
    env.update(PYTHONNOUSERSITE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    consumer = work / 'venv'
    venv.EnvBuilder(with_pip=True).create(consumer)
    python = consumer / 'bin/python'
    run([python, '-m', 'pip', 'install', f'numba-swig-bridge-rt=={args.runtime_version}',
         'numba==0.68.0', 'numpy>=2.4.1', 'scipy'], cwd=work, env=env)
    if mpi:
        mpi_wheels = list(wheels.glob('mpi4py-*.whl'))
        if len(mpi_wheels) != 1:
            parser.error('expected the source-built mpi4py wheel')
        run([python, '-m', 'pip', 'install', '--no-deps', *mpi_wheels], cwd=work, env=env)
    run([python, '-m', 'pip', 'install', '--no-deps', *pymfem_wheels], cwd=work, env=env)
    run([python, '-m', 'pip', 'check'], cwd=work, env=env)

    # Copy only fixture inputs; never copy the source mfem package.
    fixtures = work / 'fixtures'
    for directory in ('examples', 'data'):
        shutil.copytree(ROOT / directory, fixtures / directory,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.gf'))
    tests = fixtures / 'test/bridge'
    tests.mkdir(parents=True)
    for source in (ROOT / 'test/bridge').glob('test_*.py'):
        if source.name != 'test_nsb_state_build.py':
            shutil.copy2(source, tests / source.name)
    env['PATH'] = str(python.parent)
    mpi_command = []
    if mpi:
        tools = work / 'mpi-tools'
        tools.mkdir()
        original_path = os.environ.get('PATH')
        launcher = shutil.which(f'mpiexec.{args.backend}', path=original_path)
        if not launcher:
            parser.error('selected MPI launcher is unavailable')
        if args.backend == 'openmpi':
            # An absolute Open MPI launcher path can implicitly enable prefix
            # propagation, exposing /usr/bin (and SWIG) inside the workers.
            # Use the same wrapper for direct and acceptance-driver launches.
            wrapper = tools / 'mpiexec'
            wrapper.write_text(
                '#!/bin/sh\nexec ' + shlex.quote(launcher) + ' --noprefix "$@"\n'
            )
            wrapper.chmod(0o755)
        else:
            (tools / 'mpiexec').symlink_to(launcher)
        for name in ('ssh', 'rsh', 'orted', 'prted', 'hydra_pmi_proxy'):
            tool = shutil.which(name, path=original_path)
            if tool:
                (tools / name).symlink_to(tool)
        env['PATH'] += os.pathsep + str(tools)
        if args.backend == 'openmpi':
            env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'
        mpi_command = [tools / 'mpiexec', '-n', '2']

    check = work / 'check_install.py'
    check.write_text(
        "import importlib, importlib.util, pathlib, shutil, sys\n"
        "assert importlib.util.find_spec('numba_swig_bridge') is None\n"
        "assert importlib.util.find_spec('swig') is None\n"
        "assert shutil.which('swig') is None\n"
        + ("assert importlib.util.find_spec('mpi4py') is None\n" if not mpi else '')
        + "import mfem, nsb_rt\n"
        "assert pathlib.Path(mfem.__file__).resolve().is_relative_to(pathlib.Path(sys.prefix))\n"
        "assert pathlib.Path(nsb_rt.__file__).resolve().is_relative_to(pathlib.Path(sys.prefix))\n"
        f"assert nsb_rt.__version__ == {args.runtime_version!r}\n"
        "backend = importlib.import_module('mfem.' + sys.argv[1])\n"
        "registration = backend.get_bridge_registration()\n"
        "assert pathlib.Path(registration.__file__).resolve().is_relative_to(pathlib.Path(sys.prefix))\n"
        "print('Installed wheel and runtime verified:', mfem.__file__, flush=True)\n"
    )
    run([python, '-I', check, 'ser'], cwd=work, env=env)
    completed = []
    excluded = {'test_ex18p.py', 'test_ex9p.py', 'test_ex19p.py', 'test_parallel_foundation.py',
                'test_gridfunctionp.py', 'test_native_operatorsp.py'}
    for test in sorted(tests.glob('test_*.py')):
        if test.name not in excluded:
            run([python, '-I', test], cwd=work, env=env, timeout=1200)
            completed.append(test.name)
    if mpi:
        mpi_check = work / 'check_mpi.py'
        mpi_check.write_text(
            "from mpi4py import MPI\n"
            "version = MPI.Get_library_version()\n"
            f"assert {args.backend!r} in version.lower().replace(' ', ''), version\n"
            "assert MPI.COMM_WORLD.size == 2\n"
            "print('MPI rank', MPI.COMM_WORLD.rank, 'size', MPI.COMM_WORLD.size, flush=True)\n"
        )
        run([*mpi_command, python, '-I', mpi_check], cwd=work, env=env)
        run([*mpi_command, python, '-I', check, 'par'], cwd=work, env=env)
        for name in ('test_parallel_foundation.py', 'test_ex18p.py',
                     'test_gridfunctionp.py', 'test_native_operatorsp.py'):
            run([*mpi_command, python, '-I', tests / name], cwd=work, env=env)
            completed.append(name)
        # These acceptance drivers launch their own two-rank workers.
        for name in ('test_ex9p.py', 'test_ex19p.py'):
            run([python, '-I', tests / name], cwd=work, env=env, timeout=1800)
            completed.append(name)
        run([*mpi_command, python, '-I', fixtures / 'examples/bridge/ex4p.py',
             '--refine', '0', '--compare'], cwd=work, env=env)
        completed.append('ex4p.py --compare')
    (work / 'report.json').write_text(json.dumps({
        'backend': args.backend, 'runtime': args.runtime_version,
        'wheel': pymfem_wheels[0].name, 'checks': completed, 'status': 'passed',
    }, indent=2) + '\n')
    print('NSB installed-consumer validation passed:', args.backend, flush=True)


if __name__ == '__main__':
    main()
