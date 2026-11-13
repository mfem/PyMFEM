"""
Validate final wheel ownership, runtime composition, and uninstall/reinstall.

Usage: python _build_system/distribution/validate_mpich_pair.py BASE.whl ADDON.whl
       --work-dir NEW_DIRECTORY --examples-dir examples --mesh data/star.mesh
Use exact TestPyPI wheels downloaded with --no-deps --only-binary=:all: as
inputs when validating publication. Ordinary dependencies come from PyPI.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path, PurePosixPath

from select_wheel_payload import check_payload, read_wheel


def installed_paths(names, scheme, project):
    """
    Map wheel roots and wheel data schemes to their actual installation paths.
    """
    paths = set()
    for name in names:
        parts = PurePosixPath(name).parts
        if parts[0].endswith('.data'):
            category = parts[1]
            if category == 'headers':
                root = Path(scheme['data']) / 'include' / 'site' / (
                    f'python{sys.version_info.major}.{sys.version_info.minor}') / project
            else:
                root = Path(scheme[category])
            target = root.joinpath(*parts[2:])
        else:
            target = Path(scheme['platlib']) / name
        paths.add(str(target.resolve()))
    return paths


def main():
    """
    Run each binding mode in a fresh interpreter and verify file preservation.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base', type=Path)
    parser.add_argument('addon', type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--examples-dir', required=True, type=Path)
    parser.add_argument('--mesh', required=True, type=Path)
    args = parser.parse_args()
    base, addon = args.base.resolve(), args.addon.resolve()
    base_names = check_payload(base, repaired=True)
    addon_names = check_payload(addon, repaired=True, profile='mpich')
    base_info, _ = read_wheel(base)
    addon_info, _ = read_wheel(addon, profile='mpich')
    assert base_info.split('-', 1)[1] == addon_info.split('-', 1)[1], 'Wheel versions differ'
    version = base_info.split('-', 1)[1][:-len('.dist-info')]
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    constraints = work / 'candidate-wheels.txt'
    constraints.write_text(
        f'mfem @ {base.as_uri()}\nmfem-mpich @ {addon.as_uri()}\n', encoding='utf-8')
    environment = os.environ.copy()
    for key in ('PYTHONPATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH',
                'DYLD_FALLBACK_LIBRARY_PATH'):
        environment.pop(key, None)

    def run(command):
        print('Running:', ' '.join(map(str, command)), flush=True)
        return subprocess.run(list(map(str, command)), check=True, cwd=work,
                              env=environment)

    runtime = work / 'runtime'
    run([sys.executable, '-m', 'venv', runtime])
    python = runtime / 'bin' / 'python'
    pip = [python, '-m', 'pip']
    scheme = json.loads(subprocess.check_output(
        [str(python), '-I', '-c', 'import json, sysconfig; print(json.dumps(sysconfig.get_paths()))'],
        env=environment, text=True))
    overlap = installed_paths(base_names, scheme, 'mfem') & installed_paths(
        addon_names, scheme, 'mfem-mpich')
    (work / 'overlap.json').write_text(json.dumps(sorted(overlap), indent=2) + '\n')
    assert not overlap, f'Overlapping installed paths: {sorted(overlap)}'
    for name in ('check_serial_install.py', 'check_mpich_install.py'):
        shutil.copyfile(Path(__file__).with_name(name), work / name)
    for name in ('ex1.py', 'ex1p.py'):
        shutil.copyfile(args.examples_dir / name, work / name)
    shutil.copyfile(args.mesh, work / 'star.mesh')

    def smoke():
        run([runtime / 'bin/mpiexec', '-n', '2', python, '-I',
             work / 'check_mpich_install.py'])

    def snapshot(project):
        code = '''import hashlib, importlib.metadata, json
files = importlib.metadata.distribution(PROJECT).files
result = {}
for entry in files:
    path = importlib.metadata.distribution(PROJECT).locate_file(entry).resolve()
    assert path.is_file(), str(path)
    result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
print(json.dumps(result))
'''.replace('PROJECT', repr(project))
        return json.loads(subprocess.check_output(
            [str(python), '-I', '-c', code], cwd=work, env=environment, text=True))

    run(pip + ['install', '--only-binary=:all:', base])
    run([python, '-I', work / 'check_serial_install.py', work / 'ex1.py', work / 'star.mesh'])
    run(pip + ['install', '--only-binary=:all:', '--constraint', constraints,
               f'mfem[mpich] @ {base.as_uri()}'])
    run([python, '-I', '-c',
         f'from importlib.metadata import version; assert version("mfem-mpich") == {version!r}'])
    run([python, '-I', '-m', 'mfem.common.find_libmpi', '--check'])
    smoke()
    addon_before = snapshot('mfem-mpich')
    assert not set(snapshot('mfem')) & set(addon_before), 'Installed ownership overlaps'
    run(pip + ['install', '--force-reinstall', '--no-deps', base])
    assert snapshot('mfem-mpich') == addon_before, 'Base reinstall changed add-on files'
    smoke()
    run(pip + ['uninstall', '-y', 'mfem-mpich'])
    run([python, '-I', work / 'ex1.py', '-m', work / 'star.mesh'])
    run([python, '-I', '-c', "import importlib.util; assert importlib.util.find_spec('mfem.par') is None; assert importlib.util.find_spec('mfem._par') is None"])
    run(pip + ['install', '--no-deps', addon])
    addon_before = snapshot('mfem-mpich')
    run(pip + ['uninstall', '-y', 'mfem'])
    assert snapshot('mfem-mpich') == addon_before, 'Base uninstall changed add-on files'
    run(pip + ['install', '--no-deps', base])
    smoke()
    run(pip + ['check'])
    result = {name: hashlib.sha256(path.read_bytes()).hexdigest()
              for name, path in (('base_sha256', base), ('addon_sha256', addon))}
    result['status'] = 'passed'
    (work / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Wheel composition and installation lifecycle passed.', flush=True)


if __name__ == '__main__':
    main()
