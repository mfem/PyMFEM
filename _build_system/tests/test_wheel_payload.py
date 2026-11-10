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
        self.files = {name: b'# runtime\n' for name in payload.SERIAL_FILES}
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
            'mfem/external/par/lib/libmfem.so': b'parallel MFEM',
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

    def raw_wheel(self):
        wheel = self.root / 'mfem-1.0-cp312-cp312-linux_x86_64.whl'
        with WheelFile(wheel, 'w') as archive:
            for name, data in self.files.items():
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

    def test_record_corruption_is_rejected(self):
        raw = self.raw_wheel()
        corrupt = self.root / 'corrupt.whl'
        with ZipFile(raw) as source, ZipFile(corrupt, 'w') as target:
            for name in source.namelist():
                target.writestr(name, b'corrupted' if name == 'mfem/ser.py'
                               else source.read(name))
        with self.assertRaisesRegex(ValueError, 'Invalid RECORD'):
            payload.select_payload(corrupt, self.root / 'selected')


if __name__ == '__main__':
    unittest.main()
