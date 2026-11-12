"""
Validate a coordinated upgrade between two local wheel pairs.

Usage: python validate_mpich_upgrade.py PREVIOUS_WHEELS CURRENT_WHEELS NEW_WORK_DIR
Run from the repository root after both pairs have passed validate_mpich_pair.py.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from packaging.utils import parse_wheel_filename


def main():
    """
    Upgrade through the public extra, diagnose a mismatched pair, and restore it.
    """
    previous, current, work = (Path(arg).resolve() for arg in sys.argv[1:])
    old_base, = previous.glob('mfem-*.whl')
    new_base, = current.glob('mfem-*.whl')
    old_version = str(parse_wheel_filename(old_base.name)[1])
    new_version = str(parse_wheel_filename(new_base.name)[1])
    assert parse_wheel_filename(old_base.name)[1] < parse_wheel_filename(new_base.name)[1]
    work.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    for key in ('PYTHONPATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH',
                'DYLD_FALLBACK_LIBRARY_PATH'):
        environment.pop(key, None)

    def run(command):
        print('Running:', ' '.join(map(str, command)), flush=True)
        return subprocess.run(list(map(str, command)), check=True, cwd=work, env=environment)

    run([sys.executable, '-m', 'venv', work / 'runtime'])
    python = work / 'runtime/bin/python'
    pip = work / 'runtime/bin/pip'
    repo = Path(__file__).resolve().parents[2]
    shutil.copyfile(Path(__file__).with_name('check_mpich_install.py'), work / 'check_mpich_install.py')
    shutil.copyfile(repo / 'examples/ex1p.py', work / 'ex1p.py')
    shutil.copyfile(repo / 'data/star.mesh', work / 'star.mesh')

    def check(version):
        run([python, '-I', '-c',
             f'from importlib.metadata import version; '
             f'assert version("mfem") == version("mfem-mpich") == {version!r}'])
        run([python, '-I', repo / 'examples/ex1.py', '-m', work / 'star.mesh'])
        run([work / 'runtime/bin/mpiexec', '-n', '2', python, '-I', work / 'check_mpich_install.py'])
        run([pip, 'check'])

    run([pip, 'install', '--only-binary=:all:', '--find-links', previous,
         f'mfem[mpich]=={old_version}'])
    check(old_version)
    run([pip, 'install', '--upgrade', '--only-binary=:all:', '--find-links', current,
         f'mfem[mpich]=={new_version}'])
    check(new_version)
    run([pip, 'install', '--force-reinstall', '--no-deps', old_base])
    run([python, '-I', '-c', '''try:
    import mfem.par
except ImportError as error:
    assert 'Incompatible PyMFEM wheels' in str(error), str(error)
    print(error)
else:
    raise AssertionError('Mismatched pair imported successfully')
'''])
    broken = subprocess.run([str(pip), 'check'], cwd=work, env=environment)
    assert broken.returncode != 0, 'pip check accepted mismatched wheel versions'
    run([pip, 'install', '--upgrade', '--only-binary=:all:', '--find-links', current,
         f'mfem[mpich]=={new_version}'])
    check(new_version)
    (work / 'validation.json').write_text(json.dumps({
        'previous': old_version, 'current': new_version, 'status': 'passed'}, indent=2) + '\n')
    print('Coordinated upgrade and mismatch diagnostic passed.', flush=True)


if __name__ == '__main__':
    main()
