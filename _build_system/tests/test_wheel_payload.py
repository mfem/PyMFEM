"""
Release payload checks that run without compiling or importing native PyMFEM.
Run with: python -m unittest discover -s _build_system/tests
"""

import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from wheel.wheelfile import WheelFile

spec = importlib.util.spec_from_file_location(
    'select_wheel_payload',
    Path(__file__).resolve().parents[1] / 'select_wheel_payload.py')
payload = importlib.util.module_from_spec(spec)
spec.loader.exec_module(payload)


class WheelPayloadTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.dist_info = 'mfem-1.0.dist-info'
        self.files = {name: b'# runtime\n' for name in payload.PROFILES['serial']['required']}
        self.files.update({
            'mfem/_ser/__init__.py': b'',
            'mfem/_ser/mesh.py': b'# generated proxy',
            'mfem/_ser/_mesh.cpython-312-x86_64-linux-gnu.so': b'extension',
            'mfem/external/ser/lib/libmfem.so.4.10': b'native',
            'mfem/external/ser/include/mfem.hpp': b'// installed header',
            'mfem/external/ser/data/star.mesh': b'mesh',
            'mfem/par.py': b'# parallel',
            'mfem/_par/__init__.py': b'# parallel',
            'mfem/_par/_mesh.so': b'parallel extension',
            'mfem/_par/gslib.py': b'gslib proxy',
            'mfem/_par/_gslib.so': b'gslib extension',
            'mfem/external/par/lib/libmfem.so': b'parallel MFEM',
            'mfem/external/par/include/mfem.hpp': b'parallel header',
            'mfem/external/par/data/star.mesh': b'mesh',
            'mfem/external/bin/gpmetis': b'not runtime',
            'mfem/common/chypre.py': b'shared optional utility',
            'mfem/external/lib/libHYPRE.so': b'HYPRE',
            'mfem/external/lib64/libmetis.so': b'METIS',
            'mfem/_ser/setup.py': b'# build script',
            'mfem/_ser/mesh.i': b'// swig',
            'mfem/_ser/mesh_wrap.cxx': b'// generated source',
            'mfem/common/generate_bilininteg_ext.py': b'# generator',
            'mfem/common/pycoefficient.hpp': b'// wrapper source',
            self.dist_info + '/METADATA': (
                b'Metadata-Version: 2.4\nName: mfem\nVersion: 1.0\n'
                b'Requires-Dist: numpy>=2.4.1\nRequires-Dist: scipy\n'
                b'Requires-Dist: numba\nProvides-Extra: MPI\n'),
            self.dist_info + '/WHEEL': (
                b'Wheel-Version: 1.0\nRoot-Is-Purelib: false\n'
                b'Tag: cp312-cp312-linux_x86_64\n'),
            self.dist_info + '/licenses/LICENSE': b'BSD license',
        })

    def raw_wheel(self, profile='serial'):
        project = payload.PROFILES[profile]['project']
        dist_info = project.replace('-', '_') + '-1.0.dist-info'
        wheel = self.root / (project.replace('-', '_') + '-1.0-cp312-cp312-linux_x86_64.whl')
        with WheelFile(wheel, 'w') as archive:
            for name, data in self.files.items():
                if name.startswith(self.dist_info + '/'):
                    name = dist_info + name[len(self.dist_info):]
                    if name.endswith('/METADATA'):
                        data = data.replace(b'Name: mfem\n', b'Name: ' + project.encode() + b'\n')
                archive.writestr(name, data)
        return wheel

    def test_selection_preserves_runtime_metadata_and_input(self):
        raw = self.raw_wheel()
        original = hashlib.sha256(raw.read_bytes()).digest()
        selected = payload.select_payload(raw, self.root / 'selected')
        self.assertEqual(hashlib.sha256(raw.read_bytes()).digest(), original)
        names = payload.check_payload(selected)
        self.assertIn('mfem/external/ser/include/mfem.hpp', names)
        self.assertIn('mfem/common/chypre.py', names)
        self.assertNotIn('mfem/par.py', names)
        self.assertNotIn('mfem/_ser/setup.py', names)
        self.assertNotIn('mfem/common/pycoefficient.hpp', names)
        with ZipFile(selected) as archive:
            for name in names:
                if not name.endswith('/RECORD'):
                    self.assertEqual(archive.read(name), self.files[name])
        with self.assertRaises(FileExistsError):
            payload.select_payload(raw, selected.parent)

    def test_missing_required_payload_is_rejected(self):
        del self.files['mfem/external/ser/lib/libmfem.so.4.10']
        with self.assertRaisesRegex(ValueError, 'Missing serial MFEM library'):
            payload.select_payload(self.raw_wheel(), self.root / 'selected')
        self.assertFalse((self.root / 'selected').exists())

    def test_wrong_identity_is_rejected(self):
        name = self.dist_info + '/METADATA'
        self.files[name] = self.files[name].replace(b'Name: mfem\n', b'Name: mfem-mpich\n')
        with self.assertRaisesRegex(ValueError, 'project identity mfem'):
            payload.select_payload(self.raw_wheel(), self.root / 'selected')

    def test_input_directory_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'must be separate'):
            payload.select_payload(self.raw_wheel(), self.root)

    def test_check_rejects_parallel_payload_even_after_repair(self):
        with self.assertRaisesRegex(ValueError, 'outside serial ownership'):
            payload.check_payload(self.raw_wheel(), repaired=True)

    def test_repaired_private_libraries_are_allowed_only_in_check(self):
        raw = self.raw_wheel()
        selected = payload.select_payload(raw, self.root / 'selected')
        # Model repair by adding a private library and regenerating RECORD.
        repaired = self.root / 'repaired' / selected.name
        repaired.parent.mkdir()
        with ZipFile(selected) as source, WheelFile(repaired, 'w') as target:
            for name in source.namelist():
                if not name.endswith('/RECORD'):
                    target.writestr(name, source.read(name))
            target.writestr('mfem.libs/libz-123.so.1', b'private dependency')
        payload.check_payload(repaired, repaired=True)
        with self.assertRaisesRegex(ValueError, 'outside serial ownership'):
            payload.check_payload(repaired)
        with self.assertRaisesRegex(ValueError, 'before native library repair'):
            payload.select_payload(repaired, self.root / 'reselected')

    def test_both_profiles_preserve_disjoint_dynamic_runtime_files(self):
        self.files.update({
            'mfem/common/new_runtime.py': b'new shared utility',
            'mfem/new_package/__init__.py': b'new runtime package',
            'mfem/_new_extension.so': b'serial extension',
            'mfem/_par/new_proxy.py': b'new parallel proxy',
        })
        selected_names = {}
        for profile in payload.PROFILES:
            raw = self.raw_wheel(profile)
            original_hash = hashlib.sha256(raw.read_bytes()).digest()
            selected = payload.select_payload(raw, self.root / profile, profile)
            self.assertEqual(hashlib.sha256(raw.read_bytes()).digest(), original_hash)
            names = set(payload.check_payload(selected, profile=profile))
            selected_names[profile] = {n for n in names if '.dist-info/' not in n}
            with ZipFile(raw) as source, ZipFile(selected) as archive:
                for name in names:
                    if not name.endswith('/RECORD'):
                        self.assertEqual(archive.read(name), source.read(name))
            self.assertNotIn('mfem/external/bin/gpmetis', names)
        self.assertFalse(selected_names['serial'] & selected_names['mpich'])
        self.assertIn('mfem/common/new_runtime.py', selected_names['serial'])
        self.assertIn('mfem/new_package/__init__.py', selected_names['serial'])
        self.assertIn('mfem/_new_extension.so', selected_names['serial'])
        self.assertIn('mfem/_par/new_proxy.py', selected_names['mpich'])
        self.assertIn('mfem/external/lib64/libmetis.so', selected_names['mpich'])

    def test_profiles_reject_wrong_identity_and_unselected_payload(self):
        for profile, other in (('serial', 'mpich'), ('mpich', 'serial')):
            with self.subTest(profile=profile):
                with self.assertRaisesRegex(ValueError, 'project identity'):
                    payload.read_wheel(self.raw_wheel(other), profile=profile)
                with self.assertRaisesRegex(ValueError, 'outside .* ownership'):
                    payload.check_payload(self.raw_wheel(profile), profile=profile)

    def test_unknown_members_are_rejected_by_both_profiles(self):
        for name in ('unexpected.txt', 'mfem/common/unknown.dat'):
            self.files[name] = b'unclassified'
            for profile in payload.PROFILES:
                with (self.subTest(name=name, profile=profile),
                      self.assertRaisesRegex(ValueError, 'Unclassified')):
                    payload.select_payload(self.raw_wheel(profile),
                                           self.root / profile, profile)
            del self.files[name]

    def test_mpich_required_libraries_and_gslib_are_checked(self):
        for name in ('mfem/external/lib/libHYPRE.so',
                     'mfem/external/lib64/libmetis.so', 'mfem/_par/_gslib.so',
                     'mfem/external/par/lib/libmfem.so'):
            original = self.files.pop(name)
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'Missing'):
                payload.select_payload(self.raw_wheel('mpich'), self.root / 'mpich', 'mpich')
            self.files[name] = original

    def test_repair_libraries_are_profile_specific(self):
        for profile, own, other in (
                ('serial', 'mfem.libs/libz.so.1', 'mfem_mpich.libs/libz.so.1'),
                ('mpich', 'mfem_mpich.libs/libz.so.1', 'mfem.libs/libz.so.1')):
            raw = self.raw_wheel(profile)
            selected = payload.select_payload(raw, self.root / profile, profile)
            for directory, expected in ((own, True), (other, False)):
                repaired = self.root / 'repair' / selected.name
                repaired.parent.mkdir(exist_ok=True)
                with ZipFile(selected) as source, WheelFile(repaired, 'w') as target:
                    for name in source.namelist():
                        if not name.endswith('/RECORD'):
                            target.writestr(name, source.read(name))
                    target.writestr(directory, b'private dependency')
                if expected:
                    payload.check_payload(repaired, repaired=True, profile=profile)
                else:
                    with self.assertRaisesRegex(ValueError, 'outside .* ownership'):
                        payload.check_payload(repaired, repaired=True, profile=profile)
                with self.assertRaisesRegex(ValueError, 'before native library repair'):
                    payload.select_payload(repaired, self.root / 'invalid', profile)

    def test_record_corruption_is_rejected(self):
        for profile, corrupted_name in (('serial', 'mfem/ser.py'), ('mpich', 'mfem/par.py')):
            raw = self.raw_wheel(profile)
            corrupt = self.root / 'corrupt.whl'
            with ZipFile(raw) as source, ZipFile(corrupt, 'w') as target:
                for name in source.namelist():
                    target.writestr(name, b'corrupted' if name == corrupted_name
                                   else source.read(name))
            with self.subTest(profile=profile), self.assertRaisesRegex(ValueError, 'Invalid RECORD'):
                payload.select_payload(corrupt, self.root / 'selected', profile)


if __name__ == '__main__':
    unittest.main()
