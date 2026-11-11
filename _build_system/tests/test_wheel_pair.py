"""
Catch wheel paths that differ in archives but collide after installation.
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

with patch.object(sys, 'path', [str(Path(__file__).resolve().parents[1])] + sys.path):
    from validate_mpich_pair import installed_paths


class InstalledPathsTest(unittest.TestCase):
    def test_data_schemes_cannot_hide_installation_collisions(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            scheme = {
                'platlib': str(prefix / 'lib/site-packages'),
                'purelib': str(prefix / 'lib/site-packages'),
                'scripts': str(prefix / 'bin'),
                'data': str(prefix),
            }
            base = installed_paths(['mfem/common/runtime.py'], scheme, 'mfem')
            for member in ('mfem_mpich-1.data/platlib/mfem/common/runtime.py',
                           'mfem_mpich-1.data/purelib/mfem/common/runtime.py',
                           'mfem_mpich-1.data/data/lib/site-packages/mfem/common/runtime.py'):
                with self.subTest(member=member):
                    self.assertEqual(base, installed_paths([member], scheme, 'mfem-mpich'))
            self.assertEqual(
                installed_paths(['mfem-1.data/scripts/tool'], scheme, 'mfem'),
                installed_paths(['mfem_mpich-1.data/data/bin/tool'], scheme, 'mfem-mpich'))
            self.assertFalse(base & installed_paths(['mfem/_par/runtime.py'], scheme, 'mfem-mpich'))


if __name__ == '__main__':
    unittest.main()
