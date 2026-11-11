"""
Check MPI wrapper queries and parsing without a native build.
"""

import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

with patch.object(sys, 'path', [str(Path(__file__).resolve().parents[1])] + sys.path):
    import build_pymfem
    import build_utils


class MPIWrapperTest(unittest.TestCase):
    def query(self, command, responses):
        """
        Supply per-option responses and retain exact query arguments.
        """
        def run(args, **kwargs):
            self.assertNotIn('shell', kwargs)
            self.assertEqual(kwargs['timeout'], 10)
            code, stdout, stderr = responses.get(args[-1], (1, '', 'unknown option'))
            return subprocess.CompletedProcess(args, code, stdout, stderr)

        with patch.object(build_utils.subprocess, 'run', side_effect=run) as mocked:
            info = build_utils.get_mpi_wrapper_info(command)
        return info, mocked.call_args_list

    def test_openmpi_paths_quoting_order_and_other_flags(self):
        info, calls = self.query('"/opt/MPI tools/mpicxx"', {
            '--showme:compile': (0, '-I"/opt/MPI tools/include" -isystem /other -I/other', ''),
            '--showme:link': (0, ('-L/mpi/lib64 -L /mpi/lib -L/mpi/lib64 '
                                 '-lmpi_cxx -l mpi -lmpi -Wl,-rpath,/mpi/lib'), ''),
        })
        self.assertEqual(info['type'], 'openmpi')
        self.assertEqual(info['include_dirs'], ['/opt/MPI tools/include', '/other'])
        self.assertEqual(info['library_dirs'], ['/mpi/lib64', '/mpi/lib'])
        self.assertEqual(info['libraries'], ['mpi_cxx', 'mpi'])
        self.assertIn('-Wl,-rpath,/mpi/lib', info['flags'])
        self.assertEqual(calls[0].args[0], ['/opt/MPI tools/mpicxx', '--showme:compile'])

    def test_mpich_and_compiler_arguments(self):
        info, calls = self.query(['mpicxx', '-cxx=clang++'], {
            '-show': (0, 'clang++ -I/mpi/include -L/mpi/lib -lmpi', ''),
            '-v': (0, '', 'mpicxx for MPICH version 5.0.2'),
        })
        self.assertEqual(info['type'], 'mpich')
        self.assertEqual(info['include_dirs'], ['/mpi/include'])
        self.assertEqual(info['library_dirs'], ['/mpi/lib'])
        self.assertEqual(info['libraries'], ['mpi'])
        self.assertTrue(all(call.args[0][:2] == ['mpicxx', '-cxx=clang++'] for call in calls))

    def test_intel_named_and_generic_wrappers(self):
        for command, banner in [('mpiicpx', 'icpx version'),
                                ('mpicxx', 'Intel(R) MPI Library 2021')]:
            with self.subTest(command=command):
                info, _ = self.query(command, {
                    '-show': (0, 'icpx -I/intel/include -L/intel/lib -lmpi', ''),
                    '-v': (0, '', banner),
                })
                self.assertEqual(info['type'], 'intel')
                self.assertEqual(info['library_dirs'], ['/intel/lib'])

    def test_cray(self):
        info, calls = self.query('/opt/cray/bin/CC', {
            '--cray-print-opts=all': (
                0, ('-I/cray/mpi/include -L/cray/mpi/lib -lmpi '
                    '-L/cray/fabric/lib -lfabric -pthread'), ''),
        })
        self.assertEqual(info['type'], 'cray')
        self.assertEqual(info['include_dirs'], ['/cray/mpi/include'])
        self.assertEqual(info['library_dirs'], ['/cray/mpi/lib', '/cray/fabric/lib'])
        self.assertEqual(info['libraries'], ['mpi', 'fabric'])
        self.assertNotIn('-show', [call.args[0][-1] for call in calls])

    def test_unidentified_show_wrapper_is_not_assumed_mpich(self):
        info, _ = self.query('custom-mpicxx', {
            '-show': (0, 'g++ -I/custom/include -L/custom/lib -lmpi', ''),
        })
        self.assertEqual(info['type'], 'show-compatible')

    def test_incomplete_openmpi_query_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, 'Open MPI link query failed'):
            self.query('mpicxx', {'--showme:compile': (0, '-I/mpi/include', '')})

    def test_unsupported_empty_and_malformed_output(self):
        with self.assertRaisesRegex(RuntimeError, 'Unsupported MPI wrapper'):
            self.query('c++', {})
        for output, error in [('', 'No include or library flags'),
                              ('g++ -I', 'Missing value'),
                              ('g++ -L -lmpi', 'Missing value')]:
            with self.subTest(output=output), self.assertRaisesRegex(ValueError, error):
                self.query('mpicxx', {'-show': (0, output, '')})

    def test_missing_executable_and_timeout(self):
        for error in (FileNotFoundError('missing'),
                      subprocess.TimeoutExpired('mpicxx', 10)):
            with (self.subTest(error=error),
                  patch.object(build_utils.subprocess, 'run', side_effect=error),
                  self.assertRaisesRegex(RuntimeError, 'Cannot query MPI wrapper')):
                build_utils.get_mpi_wrapper_info('mpicxx')

    def test_wrapper_library_directories_reach_only_parallel_build(self):
        for serial in (True, False):
            with (self.subTest(serial=serial),
                  patch.object(build_pymfem, 'get_mpi_wrapper_info', return_value={
                      'library_dirs': ['/mpi/lib64', '/mpi/lib']}) as query,
                  patch.object(build_pymfem, 'record_mfem_sha'),
                  patch.object(build_pymfem, 'write_setup_local'),
                  patch.object(build_pymfem, 'chdir', return_value='/original'),
                  patch.object(build_pymfem.os, 'chdir'),
                  patch.object(build_pymfem.os.path, 'exists', return_value=True),
                  patch.object(build_pymfem, 'make_call') as call):
                build_pymfem.make_mfem_wrapper(serial=serial)
                command = call.call_args.args[0]
                if serial:
                    query.assert_not_called()
                    self.assertFalse(any(arg.startswith('--library-dirs=') for arg in command))
                else:
                    query.assert_called_once_with(build_pymfem.bglb.mpicxx_command)
                    self.assertIn('--library-dirs=/mpi/lib64:/mpi/lib', command)


if __name__ == '__main__':
    unittest.main()
