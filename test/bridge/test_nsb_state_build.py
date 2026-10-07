"""Source-specific state flags must not leak across parallel wrapper builds."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

from setuptools import Distribution, Extension

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '_build_system'))
from build_generatedwrapperext import Build_StateO3, state_optimization_flag


class StateBuildTests(unittest.TestCase):
    def test_supported_compilers(self):
        for banner in ('g++ (GCC) 8.5.0', 'clang version 18', 'Apple clang version 16'):
            self.assertEqual(state_optimization_flag('unix', banner), '-O3')
        for compiler, banner in (('msvc', 'Microsoft'), ('unix', 'unknown'),
                                 ('unix', 'Intel classic compiler'), ('unix', '')):
            self.assertIsNone(state_optimization_flag(compiler, banner))

    def exercise_build(self, *, fail=False, banner='g++ (GCC) 8.5.0',
                       compiler_type='unix', marked=True):
        extension = Extension('mfem._nsb_state_bindings_ext',
                              sources=['runtime/nsb_director_state.cpp', 'state_wrap.cxx'])
        if marked:
            extension.nsb_state_sources = ('runtime/nsb_director_state.cpp',)
        build = Build_StateO3(Distribution({'ext_modules': [extension]}))
        build.extensions = [extension]
        records = {}
        options = ['-O2', '-std=c++17']

        def original(obj, src, ext, cc_args, extra_postargs, pp_opts):
            records[src] = list(extra_postargs)
        build.compiler = SimpleNamespace(compiler_type=compiler_type,
                                         compiler_cxx=['c++'], _compile=original)

        def compile_extensions(self):
            # Same basename in another directory is not a runtime source.
            sources = ['runtime/nsb_director_state.cpp', 'state_wrap.cxx',
                       'other/nsb_director_state.cpp', 'ordinary_wrap.cxx']
            with ThreadPoolExecutor(max_workers=4) as workers:
                list(workers.map(lambda source: self.compiler._compile(
                    'object.o', source, '.cpp', [], options, []), sources))
            if fail:
                raise RuntimeError('compile failure')

        completed = SimpleNamespace(stdout=banner, stderr='')
        with patch('build_generatedwrapperext.subprocess.run', return_value=completed), \
             patch('build_generatedwrapperext.Build_NoDeprecationWarning.build_extensions', compile_extensions):
            if fail:
                with self.assertRaisesRegex(RuntimeError, 'compile failure'):
                    build.build_extensions()
            else:
                build.build_extensions()
        self.assertIs(build.compiler._compile, original)
        self.assertEqual(options, ['-O2', '-std=c++17'])
        return records

    def test_native_only_in_parallel_build(self):
        records = self.exercise_build()
        self.assertEqual(records.pop('runtime/nsb_director_state.cpp'),
                         ['-O2', '-std=c++17', '-O3'])
        for options in records.values():
            self.assertEqual(options, ['-O2', '-std=c++17'])

    def test_compiler_restored_after_failure(self):
        self.exercise_build(fail=True)

    def test_warning_suppression_and_native_optimization_work_together(self):
        extension = Extension('mfem._nsb_state_bindings_ext',
                              sources=['state.cpp', 'state_wrap.cxx'],
                              extra_compile_args=['-std=c++17'])
        extension.nsb_state_sources = ('state.cpp',)
        build = Build_StateO3(Distribution({'ext_modules': [extension]}))
        build.extensions = [extension]
        records = {}

        def original(obj, src, ext, cc_args, extra_postargs, pp_opts):
            records[src] = list(extra_postargs)
        build.compiler = SimpleNamespace(compiler_type='unix',
                                         compiler_cxx=['c++'], _compile=original)

        def compile_extensions(self):
            for source in extension.sources:
                self.compiler._compile('object.o', source, '.cpp', [],
                                       extension.extra_compile_args, [])

        completed = SimpleNamespace(stdout='g++ (GCC) 8.5.0', stderr='')
        with patch('build_generatedwrapperext.subprocess.run', return_value=completed), \
             patch('build_generatedwrapperext.build_ext.build_extensions', compile_extensions):
            build.build_extensions()
        for flags in records.values():
            self.assertIn('-Wno-deprecated-declarations', flags)
        self.assertEqual(records['state.cpp'][-1], '-O3')
        self.assertNotIn('-O3', records['state_wrap.cxx'])
        self.assertNotIn('-O3', extension.extra_compile_args)
        self.assertIs(build.compiler._compile, original)

    def test_fallback_and_bridge_disabled(self):
        for kwargs in ({'banner': 'unknown'}, {'compiler_type': 'msvc'}, {'marked': False}):
            with self.subTest(**kwargs):
                for options in self.exercise_build(**kwargs).values():
                    self.assertEqual(options, ['-O2', '-std=c++17'])


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]])
