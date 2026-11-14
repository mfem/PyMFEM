"""
Select the serial or MPICH release payload before native library repair.

Usage: python _build_system/distribution/select_wheel_payload.py RAW.whl --profile serial
       --output-dir selected-wheels
Check a selected or repaired artifact with --check instead of --output-dir.
After native repair, use --normalize-rpaths --output-dir FINAL_DIR.
Requires wheel; Linux uses readelf/patchelf, macOS uses delocate/Xcode tools;
does not import PyMFEM or the native build machinery.
"""

import argparse
import base64
import csv
import hashlib
import io
import json
import posixpath
import re
import shutil
import subprocess
import sys
import tempfile
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

# Build inputs are excluded separately from ownership. Installed external trees
# intentionally retain their headers, archives, configuration files and data.
BUILD_FILES = {
    'mfem/_ser/setup.py', 'mfem/_par/setup.py',
    'mfem/common/generate_bilininteg_ext.py',
    'mfem/common/generate_lininteg_ext.py',
}
BUILD_SUFFIXES = {'.i', '.h', '.hpp', '.c', '.cpp', '.cxx', '.pyc'}
MACHO_MAGIC = {b'\xcf\xfa\xed\xfe', b'\xfe\xed\xfa\xcf',
               b'\xca\xfe\xba\xbe', b'\xbe\xba\xfe\xca',
               b'\xca\xfe\xba\xbf', b'\xbf\xba\xfe\xca'}


def native_library(name):
    """
    Recognize installed shared libraries, including versioned ELF libraries.
    """
    path = PurePosixPath(name)
    return path.suffix == '.dylib' or '.so' in path.suffixes


def runtime_owner(name):
    """
    Classify dynamically discovered runtime files and explicit build inputs.
    """
    if name == 'mfem/_par/mpi-runtime.json':
        return 'mpich'
    if name.startswith('mfem/external/bin/'):
        return 'excluded'
    if name.startswith('mfem/external/ser/'):
        return 'serial'
    if name.startswith('mfem/external/'):
        return 'mpich'
    if name.startswith('mfem/'):
        if name in BUILD_FILES or PurePosixPath(name).suffix in BUILD_SUFFIXES:
            return 'excluded'
        if name.endswith('.py') or native_library(name):
            if name == 'mfem/par.py' or name.startswith('mfem/_par/'):
                return 'mpich'
            return 'serial'
    raise ValueError(f'Unclassified wheel member: {name}')


def installed_library(name, prefix, stem):
    """
    Match a shared library under an installed native library directory.
    """
    path = PurePosixPath(name)
    return (name.startswith(prefix) and path.parent.name in ('lib', 'lib64')
            and path.name.split('.')[0] == stem and native_library(name))


def required_components(mode):
    """
    Define minimum installed native components for each build profile.
    """
    prefix = 'mfem/external/' + ('ser/' if mode == 'serial' else 'par/')
    bindings = 'mfem/_ser/' if mode == 'serial' else 'mfem/_par/'
    checks = [
        (mode + ' extension modules', lambda n: n.startswith(bindings + '_')
         and n.endswith('.so')),
        (mode + ' MFEM library', lambda n: installed_library(n, prefix, 'libmfem')),
        ('installed MFEM headers', lambda n: n.startswith(prefix + 'include/')
         and n.endswith('mfem.hpp')),
        ('MFEM runtime data', lambda n: n.startswith(prefix + 'data/')),
    ]
    if mode == 'mpich':
        checks += [
            ('HYPRE library', lambda n: installed_library(n, 'mfem/external/', 'libHYPRE')),
            ('METIS library', lambda n: installed_library(n, 'mfem/external/', 'libmetis')),
            ('parallel GSLIB extension', lambda n: n.startswith('mfem/_par/_gslib.')
             and n.endswith('.so')),
        ]
    return checks


PROFILES = {
    'serial': {
        'project': 'mfem',
        'owns': lambda name: runtime_owner(name) == 'serial',
        'repair_dirs': ('mfem.libs/', 'mfem/.dylibs/'),
        'required': {'mfem/__init__.py', 'mfem/__sha__.py', 'mfem/ser.py',
                     'mfem/_ser/__init__.py', 'mfem/common/find_libmpi.py',
                     'mfem/common/check_parallel.py'},
        'components': required_components('serial'),
    },
    'mpich': {
        'project': 'mfem-mpich',
        'owns': lambda name: runtime_owner(name) == 'mpich',
        'repair_dirs': ('mfem_mpich.libs/', 'mfem_mpich.dylibs/'),
        'required': {'mfem/par.py', 'mfem/_par/__init__.py', 'mfem/_par/gslib.py',
                     'mfem/_par/mpi-runtime.json'},
        'components': required_components('mpich'),
    },
}
REPAIR_DIRS = tuple(directory for profile in PROFILES.values()
                    for directory in profile['repair_dirs'])


def selected_file(name, dist_info, profile, repaired=False):
    """
    Keep this distribution's metadata, runtime files and repaired libraries.
    """
    if name.startswith(dist_info + '/'):
        return True
    if name.startswith(REPAIR_DIRS):
        return (repaired and name.startswith(profile['repair_dirs'])
                and native_library(name))
    return profile['owns'](name)


def normalized_project(name):
    """
    Normalize Python distribution names for identity checks.
    """
    return re.sub(r'[-_.]+', '-', name).lower()


def read_wheel(wheel, profile='serial'):
    """
    Validate paths, project identity, and every RECORD entry before extraction.
    """
    definition = PROFILES[profile]
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
        directory_project, separator, version = dist_info[:-len('.dist-info')].rpartition('-')
        expected = definition['project']
        if (normalized_project(metadata.get('Name', '')) != expected or
                not separator or normalized_project(directory_project) != expected or
                metadata.get('Version') != version):
            raise ValueError(f'The {profile} profile requires project identity {expected}')
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


def elf_rpaths(path):
    """
    Read both ELF runtime-path tags, including empty search-path entries.
    """
    result = subprocess.run(['readelf', '--dynamic', '--wide', str(path)],
                            check=True, capture_output=True, text=True)
    return re.findall(r'\((RPATH|RUNPATH)\).*?\[([^\]]*)\]', result.stdout)


def origin_path(entry):
    """
    Accept origin-anchored loader paths, without accepting bare relative paths.
    """
    return any(entry == token or entry.startswith(token + '/')
               for token in ('$ORIGIN', '${ORIGIN}'))


def mpi_rpaths(name, directories, token):
    """
    Anchor package-recorded MPI directories at a native wheel member's location.
    Both paths are relative to the installed site-packages directory.
    """
    prefix = '../' * len(PurePosixPath(name).parent.parts)
    return [token + '/' + posixpath.normpath(prefix + directory) for directory in directories]


def check_wheel_rpaths(wheel, names, mpi_dirs=()):
    """
    Reject non-origin runtime paths in every ELF member, regardless of suffix.
    """
    with ZipFile(wheel) as archive, tempfile.TemporaryDirectory(prefix='pymfem-elf-') as root:
        path = Path(root) / 'member'
        for name in names:
            with archive.open(name) as stream:
                magic = stream.read(4)
                if magic != b'\x7fELF':
                    continue
                with path.open('wb') as target:
                    target.write(magic)
                    shutil.copyfileobj(stream, target)
            tags = elf_rpaths(path)
            for tag, value in tags:
                invalid = [entry for entry in value.split(':') if not origin_path(entry)]
                if invalid:
                    raise ValueError(f'Non-origin {tag} in {name}: {invalid!r}')
            values = dict(tags)
            current = values.get('RUNPATH', values.get('RPATH', '')).split(':')
            missing = set(mpi_rpaths(name, mpi_dirs, '$ORIGIN')) - set(current)
            if missing:
                raise ValueError(f'Missing MPI runtime paths in {name}: {sorted(missing)}')


def macos_paths(path):
    """
    Read Mach-O dependency names and runtime paths with native delocate tools.
    """
    from delocate.tools import get_install_names, get_rpaths

    return get_install_names(str(path)), get_rpaths(str(path))


def rm_duplicate_mfem(root):
    """
    Replace duplicate ELF linker aliases with scripts targeting their SONAME.
    Keep the versioned library for runtime loading and CMake imported targets,
    and keep the unversioned linker input for downstream builds using -lmfem.
    """
    for alias in (root / 'mfem' / 'external').rglob('libmfem.so'):
        data = alias.read_bytes()
        if not data.startswith(b'\x7fELF'):
            continue
        result = subprocess.run(['patchelf', '--print-soname', str(alias)],
                                check=True, capture_output=True, text=True)
        soname = result.stdout.strip()
        if not soname.startswith('libmfem.so.') or PurePosixPath(soname).name != soname:
            continue
        library = alias.with_name(soname)
        if library.is_file() and library.read_bytes() == data:
            alias.write_text(f'INPUT ( {soname} )\n', encoding='utf-8')
            print(f'Replaced duplicate {alias.relative_to(root)} with a linker script')


def macos_relative_path(path):
    """
    Accept loader-relative and executable-relative runtime search directories.
    """
    return path.startswith(('@loader_path/', '@executable_path/'))


def check_wheel_macos_paths(wheel, names, mpi_dirs=()):
    """
    Reject build paths in Mach-O runtime paths and dependency names.
    """
    with ZipFile(wheel) as archive, tempfile.TemporaryDirectory(prefix='pymfem-macho-') as root:
        path = Path(root) / 'member'
        for name in names:
            with archive.open(name) as stream:
                magic = stream.read(4)
                if magic not in MACHO_MAGIC:
                    continue
                path.write_bytes(magic + stream.read())
            dependencies, rpaths = macos_paths(path)
            invalid = [entry for entry in rpaths if not macos_relative_path(entry)]
            invalid += [entry for entry in dependencies if not (
                macos_relative_path(entry) or entry.startswith(
                    ('@rpath/', '/usr/lib/', '/System/Library/')))]
            if invalid:
                raise ValueError(f'Non-relative macOS loader paths in {name}: {invalid!r}')
            missing = set(mpi_rpaths(name, mpi_dirs, '@loader_path')) - set(rpaths)
            if missing:
                raise ValueError(f'Missing MPI runtime paths in {name}: {sorted(missing)}')


def check_payload(wheel, repaired=False, profile='serial', check_rpaths=True):
    """
    Reject missing components and files outside the selected ownership rules.
    """
    definition = PROFILES[profile]
    dist_info, names = read_wheel(wheel, profile)
    mpi_dirs = []
    if profile == 'mpich' and 'mfem/_par/mpi-runtime.json' in names:
        with ZipFile(wheel) as archive:
            config = json.loads(archive.read('mfem/_par/mpi-runtime.json'))
        if (config.get('schema') != 2 or not config.get('libraries') or
                not isinstance(config['libraries'], list) or
                not all(isinstance(name, str) and PurePosixPath(name).name == name
                        for name in config['libraries'])):
            raise ValueError('Invalid supplied MPI library inventory')
        mpi_dirs = config.get('library_dirs')
        if (not isinstance(mpi_dirs, list) or not mpi_dirs or
                not all(isinstance(directory, str) and directory and
                        not PurePosixPath(directory).is_absolute() for directory in mpi_dirs)):
            raise ValueError('Invalid supplied MPI library directories')
        # Strip version/type suffixes, then allow repair's optional hash suffix.
        stems = {re.split(r'\.so(?:\.|$)|(?:\.\d+)*\.dylib$', library)[0]
                 for library in config['libraries']}
        pattern = re.compile(r'(?:' + '|'.join(re.escape(stem) for stem in stems) +
                             r')(?:-[0-9a-f]+)?(?:\.|$)')
        forbidden = [name for name in names
                     if pattern.match(PurePosixPath(name).name)]
        if forbidden:
            raise ValueError('Bundled MPI libraries are forbidden:\n' + '\n'.join(forbidden))
        
    unexpected = [name for name in names
                  if not selected_file(name, dist_info, definition, repaired)]
    if unexpected:
        raise ValueError(f'Files outside {profile} ownership:\n' + '\n'.join(unexpected))
    
    missing = definition['required'] - set(names)
    if missing:
        raise ValueError(f'Missing {profile} files: ' + ', '.join(sorted(missing)))
    for description, predicate in definition['components']:
        if not any(predicate(name) for name in names):
            raise ValueError('Missing ' + description)
        
    if repaired and check_rpaths:
        check_wheel_rpaths(wheel, names, mpi_dirs)
        if sys.platform == 'darwin':
            check_wheel_macos_paths(wheel, names, mpi_dirs)
    return names


def select_payload(wheel, output_dir, profile='serial'):
    """
    Select into a temporary directory and repack with a regenerated RECORD.
    """
    wheel = Path(wheel).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir == wheel.parent:
        raise ValueError('Input wheel and output directory must be separate')
    
    definition = PROFILES[profile]
    dist_info, names = read_wheel(wheel, profile)
    if any(name.startswith(REPAIR_DIRS) for name in names):
        raise ValueError('Select payload before native library repair')
    
    with tempfile.TemporaryDirectory(prefix='pymfem-payload-') as temporary:
        root = Path(temporary)
        subprocess.run([sys.executable, '-m', 'wheel', 'unpack', str(wheel),
                        '--dest', str(root / 'unpacked')], check=True)
        unpacked, = (root / 'unpacked').iterdir()
        for name in names:
            if not selected_file(name, dist_info, definition):
                (unpacked / name).unlink()
        packed = root / 'packed'
        packed.mkdir()
        subprocess.run([sys.executable, '-m', 'wheel', 'pack', str(unpacked),
                        '--dest-dir', str(packed)], check=True)
        selected, = packed.glob('*.whl')
        check_payload(selected, profile=profile)
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / selected.name
        # Never overwrite the input or an earlier release artifact.
        with destination.open('xb') as stream, selected.open('rb') as source:
            shutil.copyfileobj(source, stream)
        return destination


def normalize_rpaths(wheel, output_dir, profile='serial'):
    """
    Replace build paths with relative internal and pip MPI paths after repair.
    """
    if sys.platform not in ('linux', 'darwin'):
        raise ValueError('Runtime-path normalization requires Linux or macOS')
    
    wheel, output_dir = Path(wheel).resolve(), Path(output_dir).resolve()
    if output_dir == wheel.parent:
        raise ValueError('Input wheel and output directory must be separate')
    
    names = check_payload(wheel, repaired=True, profile=profile, check_rpaths=False)
    mpi_dirs = []
    if profile == 'mpich':
        with ZipFile(wheel) as archive:
            mpi_dirs = json.loads(archive.read('mfem/_par/mpi-runtime.json'))['library_dirs']

    with tempfile.TemporaryDirectory(prefix='pymfem-rpaths-') as temporary:
        root = Path(temporary)
        subprocess.run([sys.executable, '-m', 'wheel', 'unpack', str(wheel),
                        '--dest', str(root / 'unpacked')], check=True)
        unpacked, = (root / 'unpacked').iterdir()
        changed = 0
        for name in names:
            path = unpacked / name
            with path.open('rb') as stream:
                magic = stream.read(4)
            if magic in MACHO_MAGIC and sys.platform == 'darwin':
                _, rpaths = macos_paths(path)
                invalid = set(entry for entry in rpaths if not macos_relative_path(entry))
                added = [entry for entry in mpi_rpaths(name, mpi_dirs, '@loader_path')
                         if entry not in rpaths]
                for entry in invalid:
                    subprocess.run(['install_name_tool', '-delete_rpath', entry, str(path)], check=True)
                for entry in added:
                    subprocess.run(['install_name_tool', '-add_rpath', entry, str(path)], check=True)
                if invalid or added:
                    subprocess.run(['codesign', '--force', '--sign', '-', str(path)], check=True)
                    changed += 1
                continue
            if magic != b'\x7fELF':
                continue
            tags = elf_rpaths(path)
            # RUNPATH supersedes RPATH when both tags are present. Retain the
            # effective search order rather than activating an ignored RPATH.
            values = dict(tags)
            current = values.get('RUNPATH', values.get('RPATH', ''))
            kept = [entry for entry in current.split(':') if origin_path(entry)]
            kept += [entry for entry in mpi_rpaths(name, mpi_dirs, '$ORIGIN') if entry not in kept]
            kept = ':'.join(kept)
            if kept == current and all(origin_path(entry) for _, value in tags
                                       for entry in value.split(':')):
                continue
            subprocess.run(['patchelf', '--remove-rpath', str(path)], check=True)
            if kept:
                subprocess.run(['patchelf', '--set-rpath', kept, str(path)], check=True)
            changed += 1

        # linux wheel contains libmfem and libmfem.so.(with version number)
        # in order to keep the wheel size below 100MB, we remove one.
        if sys.platform == 'linux':
            rm_duplicate_mfem(unpacked)

        packed = root / 'packed'
        packed.mkdir()
        subprocess.run([sys.executable, '-m', 'wheel', 'pack', str(unpacked),
                        '--dest-dir', str(packed)], check=True)
        normalized, = packed.glob('*.whl')
        check_payload(normalized, repaired=True, profile=profile)
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / normalized.name
        with destination.open('xb') as target, normalized.open('rb') as source:
            shutil.copyfileobj(source, target)
        print(f'Normalized runtime paths in {changed} native members')
        return destination


def main():
    """
    Select a raw wheel or check the final repaired wheel's ownership and RECORD.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('--profile', required=True, choices=PROFILES)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output-dir', type=Path)
    action.add_argument('--check', action='store_true')
    parser.add_argument('--normalize-rpaths', action='store_true',
                        help='With --output-dir, normalize an already repaired native wheel')
    args = parser.parse_args()
    if args.normalize_rpaths and args.check:
        parser.error('--normalize-rpaths requires --output-dir')
    try:
        if args.check:
            names = check_payload(args.wheel, repaired=True, profile=args.profile)
            print(f'Validated {args.profile} payload: {len(names)} files in {args.wheel}')
        else:
            operation = normalize_rpaths if args.normalize_rpaths else select_payload
            print(operation(args.wheel, args.output_dir, profile=args.profile))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
