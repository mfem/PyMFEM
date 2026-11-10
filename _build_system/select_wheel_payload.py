"""
Select the serial release payload before auditwheel/delocate repair.

Usage: python _build_system/select_wheel_payload.py RAW.whl --profile serial
       --output-dir selected-wheels
Check a selected or repaired artifact with --check instead of --output-dir.
Requires wheel; does not import PyMFEM or the native build machinery.
"""

import argparse
import base64
import csv
import hashlib
import io
import shutil
import subprocess
import sys
import tempfile
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

# Shared runtime utilities remain base-owned, including optional parallel paths.
SERIAL_FILES = {
    'mfem/__init__.py',
    'mfem/__sha__.py',
    'mfem/ser.py',
    'mfem/common/__init__.py',
    'mfem/common/arg_parser.py',
    'mfem/common/bessel.py',
    'mfem/common/chypre.py',
    'mfem/common/findpoints.py',
    'mfem/common/matrix_file.py',
    'mfem/common/mpi_debug.py',
    'mfem/common/mpi_dtype.py',
    'mfem/common/numba_coefficient_utils.py',
    'mfem/common/parcsr_extra.py',
    'mfem/common/sparse_utils.py',
}


def serial_file(name, dist_info, repaired=False):
    """
    Define base ownership; installed serial headers and data are intentional.
    """
    path = PurePosixPath(name)
    if name.startswith(dist_info + '/'):
        return True
    if name in SERIAL_FILES:
        return True
    if path.parent == PurePosixPath('mfem/_ser'):
        return path.name != 'setup.py' and (
            path.suffix == '.py' or path.name.endswith('.so'))
    if name.startswith('mfem/external/ser/'):
        return True
    # Repair tools own their relocated native dependencies. Neither directory
    # is accepted from a raw wheel, which must be selected before repair.
    if repaired and name.startswith(('mfem.libs/', 'mfem/.dylibs/')):
        return path.suffix == '.dylib' or '.so' in path.suffixes
    return False


def read_wheel(wheel):
    """
    Validate paths, project identity, and every RECORD entry before extraction.
    """
    with ZipFile(wheel) as archive:
        names = [item.filename for item in archive.infolist()
                 if not item.is_dir()]
        if len(names) != len(set(names)):
            raise ValueError('Duplicate wheel members')
        for name in names:
            path = PurePosixPath(name)
            if (path.is_absolute() or '..' in path.parts or '\\' in name or
                    path.as_posix() != name):
                raise ValueError(f'Unsafe wheel member: {name}')
        roots = {PurePosixPath(name).parts[0] for name in names
                 if PurePosixPath(name).parts[0].endswith('.dist-info')}
        if len(roots) != 1:
            raise ValueError('Expected exactly one .dist-info directory')
        dist_info = roots.pop()
        metadata = BytesParser().parsebytes(archive.read(dist_info + '/METADATA'))
        if metadata['Name'] != 'mfem' or not dist_info.startswith('mfem-'):
            raise ValueError('The serial profile requires project identity mfem')
        for required in ('METADATA', 'WHEEL', 'RECORD'):
            if dist_info + '/' + required not in names:
                raise ValueError(f'Missing {required}')
        if any(name.endswith(('/RECORD.jws', '/RECORD.p7s')) for name in names):
            raise ValueError('Select payload before signing the wheel')
        record = dist_info + '/RECORD'
        rows = list(csv.reader(io.StringIO(archive.read(record).decode('utf-8'))))
        if len(rows) != len(names) or {row[0] for row in rows} != set(names):
            raise ValueError('RECORD does not match wheel members')
        for name, digest, size in rows:
            if name == record:
                if digest or size:
                    raise ValueError('RECORD must have an empty self-entry')
                continue
            data = archive.read(name)
            algorithm, encoded = digest.split('=', 1)
            if algorithm not in ('sha256', 'sha384', 'sha512'):
                raise ValueError(f'Unsupported RECORD hash: {algorithm}')
            actual = base64.urlsafe_b64encode(
                hashlib.new(algorithm, data).digest()).rstrip(b'=').decode('ascii')
            if encoded != actual or str(len(data)) != size:
                raise ValueError(f'Invalid RECORD hash/size: {name}')
        return dist_info, names


def check_payload(wheel, repaired=False):
    """
    Reject missing serial components and files outside the ownership list.
    """
    dist_info, names = read_wheel(wheel)
    unexpected = [name for name in names
                  if not serial_file(name, dist_info, repaired)]
    if unexpected:
        raise ValueError('Files outside serial ownership:\n' + '\n'.join(unexpected))
    required = SERIAL_FILES | {'mfem/_ser/__init__.py'}
    missing = required - set(names)
    if missing:
        raise ValueError('Missing serial files: ' + ', '.join(sorted(missing)))
    for description, predicate in (
        ('serial extension modules', lambda n: n.startswith('mfem/_ser/_') and
         n.endswith('.so')),
        ('serial MFEM library', lambda n: n.startswith('mfem/external/ser/lib') and
         (PurePosixPath(n).name.startswith('libmfem') and
          ('.so' in PurePosixPath(n).suffixes or n.endswith('.dylib')))),
        ('installed MFEM headers', lambda n: n.startswith('mfem/external/ser/include/')
         and n.endswith('mfem.hpp')),
        ('MFEM runtime data', lambda n: n.startswith('mfem/external/ser/data/')),
    ):
        if not any(predicate(name) for name in names):
            raise ValueError('Missing ' + description)
    return names


def select_payload(wheel, output_dir):
    """
    Select into a temporary directory and repack with a regenerated RECORD.
    """
    wheel = Path(wheel).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir == wheel.parent:
        raise ValueError('Input wheel and output directory must be separate')
    dist_info, names = read_wheel(wheel)
    if any(name.startswith(('mfem.libs/', 'mfem/.dylibs/')) for name in names):
        raise ValueError('Select payload before native library repair')
    with tempfile.TemporaryDirectory(prefix='pymfem-payload-') as temporary:
        root = Path(temporary)
        subprocess.run([sys.executable, '-m', 'wheel', 'unpack', str(wheel),
                        '--dest', str(root / 'unpacked')], check=True)
        unpacked, = (root / 'unpacked').iterdir()
        for name in names:
            if not serial_file(name, dist_info):
                (unpacked / name).unlink()
        packed = root / 'packed'
        packed.mkdir()
        subprocess.run([sys.executable, '-m', 'wheel', 'pack', str(unpacked),
                        '--dest-dir', str(packed)], check=True)
        selected, = packed.glob('*.whl')
        check_payload(selected)
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / selected.name
        # Never overwrite the input or an earlier release artifact.
        with destination.open('xb') as stream, selected.open('rb') as source:
            shutil.copyfileobj(source, stream)
        return destination


def main():
    """
    Select a raw wheel or check the final repaired wheel's ownership and RECORD.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('--profile', required=True, choices=['serial'])
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output-dir', type=Path)
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    try:
        if args.check:
            names = check_payload(args.wheel, repaired=True)
            print(f'Validated serial payload: {len(names)} files in {args.wheel}')
        else:
            print(select_payload(args.wheel, args.output_dir))
    except (ValueError, OSError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
