"""
Check matching wheel coverage before publication.

Usage: python check_release.py DIRECTORY [--platform macosx_14_0_arm64]
"""

import argparse
from email.parser import BytesParser
from pathlib import Path
from zipfile import ZipFile

from packaging.requirements import Requirement
from packaging.utils import parse_wheel_filename


def check_release(directory, platform='manylinux_2_28_x86_64'):
    """
    Require a single release with one base/add-on pair per supported Python.
    """
    directory = Path(directory)
    coverage = {'mfem': set(), 'mfem-mpich': set()}
    versions = set()
    for wheel in directory.glob('*.whl'):
        project, version, _, tags = parse_wheel_filename(wheel.name)
        assert project in coverage, wheel.name
        versions.add(str(version))
        with ZipFile(wheel) as archive:
            member, = (name for name in archive.namelist() if name.endswith('.dist-info/METADATA'))
            metadata = BytesParser().parsebytes(archive.read(member))
        assert metadata['Version'] == str(version), wheel.name
        requirements = [Requirement(value) for value in metadata.get_all('Requires-Dist', [])]
        if project == 'mfem':
            assert 'mpich' in metadata.get_all('Provides-Extra', []), wheel.name
            peer = [r for r in requirements if r.name == 'mfem-mpich'
                    and r.marker and r.marker.evaluate({'extra': 'mpich'})
                    and not r.marker.evaluate({'extra': ''})]
        else:
            peer = [r for r in requirements if r.name == 'mfem' and r.marker is None]
        assert len(peer) == 1 and str(peer[0].specifier) == f'=={version}', wheel.name
        assert tags.isdisjoint(coverage[project]), f'Duplicate wheel coverage: {wheel.name}'
        coverage[project].update(tags)
    assert len(versions) == 1, f'Expected one release version: {versions}'
    assert coverage['mfem'] == coverage['mfem-mpich'], 'Base/add-on wheel coverage differs'
    tags = coverage['mfem']
    assert {tag.interpreter for tag in tags} == {'cp311', 'cp312', 'cp313', 'cp314'}, tags
    assert all(tag.abi == tag.interpreter for tag in tags), tags
    if platform.startswith('macosx_'):
        assert all(tag.platform == platform for tag in tags), tags
    else:
        assert all(tag.platform.startswith('manylinux_') and tag.platform.endswith('_x86_64')
                   for tag in tags), tags
    assert {tag.interpreter for tag in tags if tag.platform == platform} == {
        'cp311', 'cp312', 'cp313', 'cp314'}, f'Missing {platform} coverage'
    version, = versions
    expected_sdists = set() if platform.startswith('macosx_') else {f'mfem-{version}.tar.gz'}
    assert {path.name for path in directory.glob('*.tar.gz')} == expected_sdists
    assert not list(directory.glob('*.zip')), 'Only the base tar.gz sdist is published'
    print(f'Release {version}: matching {platform} wheel pairs passed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--platform', default='manylinux_2_28_x86_64')
    args = parser.parse_args()
    check_release(args.directory, args.platform)
